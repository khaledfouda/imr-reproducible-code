#' MCCI: Matrix Completion with Covariates and Incomplete data

source("./code/other_models/MCCI_original.R")

#' Fit MCCI with fixed hyperparameters (modified from original SMCfit)
#'
#' Uses uniform theta estimation and no survey weights (diagD = I).
#' The core SVT step uses SVTE_alpha from the original code.
#'
#' @param Aobs Observed matrix with 0s for missing entries
#' @param X Covariate matrix (n1 x m)
#' @param lambda_1 Ridge penalty for beta (tau_beta_ratio in original)
#' @param lambda_2 SVT penalty for B (tau_svd_ratio in original)
#' @param alpha Elastic net mixing parameter (alpha_ratio in original)
#' @param n1n2_optimized If TRUE, use leading eigenvalue as scaling factor
#' @param return_rank If TRUE, compute and return the rank of A_hat
#' @return List with estimates, M (B_hat), beta, and rank
MCCI_fit <- function(Aobs, X, lambda_1, lambda_2, alpha,
                     n1n2_optimized = TRUE,
                     return_rank = TRUE) {

  Aobs[is.na(Aobs)] <- 0
  n1 <- dim(Aobs)[1]
  n2 <- dim(Aobs)[2]
  m  <- dim(X)[2]

  # Observation mask
  omega <- matrix(as.numeric(Aobs != 0), n1, n2)

  # --- Projection matrices (original SMCfit L224-227 with diagD = I) ---
  X.X <- t(X) %*% X
  P_X <- X %*% MASS::ginv(X.X) %*% t(X)
  P_bar_X <- diag(1, n1) - P_X

  # --- Scaling factor (original SMCfit_cv uses svd(t(Xadj)%*%Xadj)$d[1]) ---
  if (n1n2_optimized) {
    n1n2 <- svd(X.X)$d[1]
  } else {
    n1n2 <- n1 * n2 / 2
  }

  # --- Uniform theta estimation (replaces original's logistic GLM) ---
  # Under MCAR: theta_hat = P(observed) = nnz / (n1*n2), applied uniformly.
  # We multiply by 1/theta to get inverse propensity weights.
  theta_inv <- (n1 * n2) / sum(omega)
  W_theta_Y <- Aobs * theta_inv

  # --- Beta estimation (original SMCfit L241-242, with Xadj = X since diagD = I) ---
  beta_hat <- MASS::ginv(X.X + n1n2 * lambda_1 * diag(1, m)) %*% t(X) %*% W_theta_Y

  # --- SVD of projected weighted observations (original SMCfit L238, L245) ---
  PXpAobsw <- P_bar_X %*% W_theta_Y
  svdd <- svd(PXpAobsw)

  # --- Update scaling factor for B_hat (same logic as original) ---
  if (n1n2_optimized) {
    n1n2_svd <- svdd$d[1]
  } else {
    n1n2_svd <- n1 * n2 / 2
  }

  # --- SVT for B_hat (uses original's SVTE_alpha, L247-250 / L27-32) ---
  # SVTE_alpha(u, d, v, tau_svd_ratio, alpha_ratio) computes:
  #   u %*% (pmax(d - d[1]*tau*alpha, 0) * t(v)) / (1 + 2*d[1]*tau*(1-alpha))
  B_hat <- SVTE_alpha(svdd$u, svdd$d, svdd$v,
                      tau_svd_ratio = lambda_2,
                      alpha_ratio = alpha)

  # --- Final estimate (original SMCfit L258, without sqrt(diagD) since diagD = I) ---
  A_hat <- X %*% beta_hat + B_hat

  rank <- NULL
  if (return_rank)
    rank <- qr(A_hat)$rank

  return(list(
    estimates = A_hat,
    M = B_hat,
    beta = beta_hat,
    rank = rank
  ))
}


#' Optimized fit for cross-validation (avoids recomputing projections)
#'
#' Pre-computed data from prepare_fold_data is reused across grid points.
MCCI_fit_optimized <- function(data, lambda_1, lambda_2, alpha) {
  beta_hat <- MASS::ginv(data$X.X + data$n1n2Im * lambda_1) %*% data$X.W.theta.Y
  B_hat <- SVTE_alpha(data$svdd$u, data$svdd$d, data$svdd$v,
                      tau_svd_ratio = lambda_2,
                      alpha_ratio = alpha)
  A_hat <- data$X %*% beta_hat + B_hat
  return(A_hat[data$W_fold == 0 & data$W == 1])
}


#' Pre-compute fold-specific data for efficient CV
prepare_fold_data <- function(Y_train, Y_valid, W_fold, W, X, n1n2_optimized) {
  n1 <- dim(Y_train)[1]
  n2 <- dim(Y_train)[2]
  m  <- dim(X)[2]

  X.X <- t(X) %*% X
  P_X <- X %*% MASS::ginv(X.X) %*% t(X)
  P_bar_X <- diag(1, n1) - P_X

  # Uniform theta on training fold
  theta_inv <- (n1 * n2) / sum(W_fold)
  W_theta_Y <- Y_train * theta_inv

  if (n1n2_optimized) {
    n1n2Im <- svd(X.X)$d[1] * diag(1, m)
  } else {
    n1n2Im <- n1 * n2 * diag(1, m) * 0.5
  }

  X.W.theta.Y <- t(X) %*% W_theta_Y
  svdd <- svd(P_bar_X %*% W_theta_Y)

  list(
    Y_valid = Y_valid,
    W_fold = W_fold,
    W = W,
    X = X,
    X.X = X.X,
    n1n2Im = n1n2Im,
    X.W.theta.Y = X.W.theta.Y,
    svdd = svdd
  )
}


#' K-fold cell-level split for matrix cross-validation
#'
#' Splits observed cells into folds, balanced per column.
MC_Kfold_split <- function(n_rows, n_cols, n_folds, obs_mask, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)
  indices <- indices[obs_mask == 1, ]
  indices <- indices[sample(1:nrow(indices)), ]

  indices <- indices %>%
    mutate(row = as.numeric(row), col = as.numeric(col)) %>%
    group_by(col) %>%
    do(sample_n(., size = nrow(.))) %>%
    mutate(fold = rep(1:n_folds, length.out = n())) %>%
    ungroup()

  folds <- vector("list", n_folds)
  for (i in 1:n_folds) {
    valid_mask <- matrix(1, nrow = n_rows, ncol = n_cols)
    test_indices <- indices[indices$fold == i, ]
    valid_mask[obs_mask == 0] <- 1
    valid_mask[as.matrix(test_indices[, c("row", "col")])] <- 0
    folds[[i]] <- valid_mask
  }
  return(folds)
}


#' Cross-validation for MCCI hyperparameters
#'
#' Searches over (lambda_1, lambda_2, alpha) using parallel grid search.
#' Default grid values match the original SMCfit_cv implementation.
#'
#' @param Y Response matrix (missing entries as NA or 0)
#' @param X Covariate matrix
#' @param W Binary observation mask (1 = observed)
#' @param n_folds Number of CV folds (default 5, matching original)
#' @param lambda_1_grid Grid for ridge penalty (default seq(0,1,length=30), matching original tau1_grid)
#' @param lambda_2_grid Grid for SVT penalty (default seq(0.9,0.1,length=30), matching original tau2_grid)
#' @param alpha_grid Grid for elastic net mixing (default seq(0.992,1,length=20), matching original)
#' @param seed Random seed
#' @param numCores Number of parallel cores
#' @param n1n2_optimized Use eigenvalue-based scaling factor
#' @param test_error Error metric function(predicted, true)
#' @param return_diagn If TRUE, return full grid results
#' @return List with best_parameters, best_score, fit, timing info
MCCI.cv <- function(Y, X, W,
                    n_folds = 5,
                    lambda_1_grid = seq(0, 1, length = 30),
                    lambda_2_grid = seq(0.9, 0.1, length = 30),
                    alpha_grid = seq(0.992, 1, length = 20),
                    seed = NULL,
                    numCores = 1,
                    n1n2_optimized = TRUE,
                    test_error = IMR::get_metric("rmse"),
                    return_diagn = FALSE) {

  start_time <- Sys.time()
  fit_counter <- 0

  Y[is.na(Y)] <- 0
  if (!is.null(seed)) set.seed(seed)

  best_score <- Inf
  best_params <- list(alpha = NA, lambda_1 = NA, lambda_2 = NA)

  folds <- MC_Kfold_split(nrow(Y), ncol(Y), n_folds, W, seed)

  fold_data <- lapply(1:n_folds, function(i) {
    W_fold <- folds[[i]]
    Y_train <- Y
    Y_train[W_fold == 0] <- 0
    Y_valid <- Y[W_fold == 0 & W == 1]
    prepare_fold_data(Y_train, Y_valid, W_fold, W, X, n1n2_optimized)
  })

  # Parallel grid search
  require(parallel)
  require(doParallel)
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("MCCI_fit_optimized", "SVTE_alpha"))

  results <-
    foreach(alpha = alpha_grid, .combine = rbind) %:%
    foreach(lambda_1 = lambda_1_grid, .combine = rbind) %:%
    foreach(lambda_2 = lambda_2_grid, .combine = rbind) %dopar% {
      score <- 0
      for (i in 1:n_folds) {
        data <- fold_data[[i]]
        A_hat_test <- MCCI_fit_optimized(data, lambda_1, lambda_2, alpha)
        score <- score + test_error(A_hat_test, data$Y_valid)
        fit_counter <- fit_counter + 1
      }
      score <- score / n_folds
      c(alpha, lambda_1, lambda_2, score)
    }

  # Select best: minimum score, tie-break by highest lambda_2
  min_score <- min(results[, 4])
  min_results <- results[results[, 4] == min_score, , drop = FALSE]
  if (nrow(min_results) > 1) {
    best_result <- min_results[which.max(min_results[, 3]), ]
  } else {
    best_result <- min_results
  }
  best_params <- list(
    alpha = best_result[1],
    lambda_1 = best_result[2],
    lambda_2 = best_result[3]
  )
  best_score <- best_result[4]
  stopCluster(cl)

  # Final refit on full data
  best_fit <- MCCI_fit(
    Aobs = Y, X = X,
    lambda_1 = best_params$lambda_1,
    lambda_2 = best_params$lambda_2,
    alpha = best_params$alpha,
    n1n2_optimized = n1n2_optimized,
    return_rank = TRUE
  )
  fit_counter <- fit_counter + 1

  obj <- list(
    best_parameters = best_params,
    best_score = best_score,
    fit = best_fit,
    total_num_fits = fit_counter,
    time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  )
  obj$time_per_fit <- obj$time / obj$total_num_fits
  if (return_diagn) obj$results <- results
  return(obj)
}
