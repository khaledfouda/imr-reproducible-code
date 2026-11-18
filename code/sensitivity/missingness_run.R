# Sensitivity analysis, missingness mechanism: run the simulations.
#
# One replication:
#   1. a complete data set with n_max rows, theta = X beta + Gamma Z' + A B'
#      and Y = theta + noise, where conditions (C1)-(C3) hold on the n_max
#      rows;
#   2. seven masks on it, all drawn from the same uniforms: MCAR, three MAR
#      mechanisms of increasing strength driven by X, and three driven by X
#      and Z equally;
#   3. for each mask and each n, a fit on the first n rows, with lambda_m
#      tuned and lambda_beta = lambda_gamma = 0.
# The fits of a replication share everything except the mask and n. The
# errors of the fitted matrix and of Gamma are always measured on the first
# min(n) rows.
#
# Run from the root of the simulation repository, in a terminal, with one
# BLAS thread per worker:
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/missingness_run.R lite
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/missingness_run.R draft
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/missingness_run.R full
# A run that is stopped and started again continues where it stopped: each
# replication has its own result file. Tables and figure:
# missingness_results.R.

size <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(size)) size <- "lite"
stopifnot(size %in% c("lite", "draft", "full"))

# Setup ----
library(IMR)
source("./code/simulation/simulation_function.R")
# run_tasks() and collect_tasks(): parallel workers, one result file per task
source("./code/sensitivity/runner.R")



# Settings ----
cfg <- list(
  n = c(200, 400, 600, 800, 1000),  # numbers of rows; the fits are nested
  m = 400, p = 4, q = 4, r = 4, snr = 3,   # q = p
  rate = 0.3,                   # mean probability that an entry is observed
  # The mechanisms: logit(prob_ij) = a0 + kappa_x * score_i + kappa_z * score_j,
  # where score_i comes from the X of row i and score_j from the Z of column
  # j. Both zero: MCAR. kappa_z = 0: MAR driven by X. kappa_z = kappa_x: MAR
  # driven by X and Z equally.
  kappa_x = c(0, 0.5, 0.75, 1, 0.25, 0.5, 0.75),
  kappa_z = c(0, 0,   0,    0, 0.25, 0.5, 0.75),
  val_prop = 0.2, thresh = 1e-6, maxit = 2000,
  reps = 500, n_workers = 8,
  out_dir = "./output/sensitivity"
)
if (size == "draft") cfg$reps <- 27
if (size == "lite") {
  cfg$n <- c(100, 200, 300)
  cfg$m <- 200
  cfg$reps <- 2
  cfg$n_workers <- 2
  cfg$out_dir <- "./output/sensitivity_lite"
}
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  warning("OPENBLAS_NUM_THREADS is not 1: each worker may use every core.")
}

# Relative error ||a - b||_F / ||b||_F
rel_err <- function(a, b) {
  return(sqrt(sum((a - b)^2) / sum(b^2)))
}

# Data ----
# The complete data set of a replication: X, Z, beta, Gamma, theta and Y
# with n_max rows, the uniforms behind every mask, the MAR score of each
# row, and the probability that each entry is observed under every
# mechanism (a list of n_max x m matrices).
make_data <- function(seed, cfg) {
  n_max <- max(cfg$n)
  # centered: each column of X and of Z has mean zero (X on the n_max rows).
  # identified: conditions (C1)-(C3), X'A = 0, Z'B = 0 and X'Gamma = 0, hold
  # on the n_max rows. On the first n < n_max rows, the centering of X and
  # the two conditions that involve X hold only approximately.
  d <- generate_simulated_data(
    n = n_max, m = cfg$m, r = cfg$r, p = cfg$p, q = cfg$q, seed = seed,
    centered = TRUE, identified = TRUE
  )
  theta <- d$theta   # X beta + Gamma Z' + A B'

  # noise, direction of the row score and uniforms, with a seed of their own.
  # The direction of the column score is drawn last, so that the draws
  # before it do not depend on it.
  set.seed(1e6 + seed)
  sigma <- sd(theta) / cfg$snr
  noise <- matrix(rnorm(n_max * cfg$m, sd = sigma), n_max, cfg$m)
  w <- rnorm(cfg$p)
  u <- matrix(runif(n_max * cfg$m), n_max, cfg$m)
  w_z <- rnorm(cfg$q)

  # The scores are the standardized X and Z times random unit vectors. For
  # each mechanism, a0 gives mean(prob) = rate over the n_max x m entries.
  score <- drop(scale(d$X) %*% (w / sqrt(sum(w^2))))
  score_z <- drop(scale(d$Z) %*% (w_z / sqrt(sum(w_z^2))))
  prob <- lapply(seq_along(cfg$kappa_x), function(k) {
    logit <- outer(cfg$kappa_x[k] * score, cfg$kappa_z[k] * score_z, "+")
    a0 <- uniroot(
      function(a) mean(plogis(a + logit)) - cfg$rate,
      interval = c(-10, 10), tol = 1e-10
    )$root
    return(plogis(a0 + logit))
  })

  return(list(
    X = d$X, Z = d$Z, beta = d$beta, gamma = d$gamma, theta = theta,
    Y = theta + noise, u = u, prob = prob, score = score
  ))
}

# One fit ----
# The first n rows under mechanism k: the tuning run and the errors, as one
# row of results.
fit_cell <- function(d, n, k, seed, cfg) {
  rows <- seq_len(n)
  test <- seq_len(min(cfg$n))   # the rows that every size contains
  x_mat <- d$X[rows, ]
  # the first n rows of the mask of the complete data set (TRUE = observed),
  # so the masks of the different sizes are nested
  prob <- d$prob[[k]][rows, ]
  mask <- d$u[rows, ] < prob
  # the identified truth of these rows, T(theta): beta = (X'X)^{-1} X' theta
  # and Gamma = Px_perp theta Z (Z'Z)^{-1}. They are the generated beta and
  # Gamma when X'A = 0 and X'Gamma = 0 on these rows, as at n = n_max.
  qx <- qr(x_mat)
  beta_t <- qr.coef(qx, d$theta[rows, ])
  gamma_t <- t(qr.coef(qr(d$Z), t(qr.resid(qx, d$theta[rows, ]))))
  cell <- data.frame(
    seed = seed, n = n, kappa_x = cfg$kappa_x[k], kappa_z = cfg$kappa_z[k],
    pi_ratio = max(prob) / min(prob),
    obs_rate = mean(mask),
    min_row_obs = min(rowSums(mask)), min_col_obs = min(colSums(mask)),
    # their distances to the generated beta and Gamma (Gamma: test rows)
    beta_gap = rel_err(beta_t, d$beta),
    gamma_gap = rel_err(gamma_t[test, ], d$gamma[test, ])
  )
  # imr() needs two observed entries in every row and column
  if (min(cell$min_row_obs, cell$min_col_obs) < 2) {
    return(cbind(cell, status = "skipped"))
  }

  # lambda_m is tuned on a validation set. The rank is fixed, and beta and
  # Gamma are dense: lambda_beta = lambda_gamma = 0.
  start <- Sys.time()
  out <- suppressMessages(imr(
    d$Y[rows, ] * mask, x_mat, d$Z,
    center = FALSE, scale = TRUE, val_prop = cfg$val_prop,
    grid = imr_tune_grid(beta = 0, gamma = 0, rank = cfg$r),
    maxit = cfg$maxit, thresh = cfg$thresh,
    n_cores = 1, seed = seed, verbose = 0
  ))
  secs <- as.numeric(Sys.time() - start, units = "secs")

  # raw estimates, and identified ones: T(theta hat)
  raw <- reconstruct(out$fit, out$data, trace = FALSE)
  id <- reconstruct(imr_identify(out$fit, out$data), out$data, trace = FALSE)

  err2 <- (raw$estimates[test, ] - d$theta[test, ])^2
  theta2 <- d$theta[test, ]^2
  unobs <- !mask[test, ]
  # the tenths of the test rows with the smallest and the largest score: the
  # rows that the MAR mechanisms observe least and most
  score_rank <- rank(d$score[test])
  tenth <- length(test) / 10
  low <- which(score_rank <= tenth)
  high <- which(score_rank > length(test) - tenth)
  return(cbind(
    cell, status = "ok",
    # fitted matrix on the test rows: all entries, the unobserved ones, and
    # all entries of the least and of the most observed rows
    theta_all = sqrt(sum(err2) / sum(theta2)),
    theta_unobs = sqrt(sum(err2[unobs]) / sum(theta2[unobs])),
    theta_low = sqrt(sum(err2[low, ]) / sum(theta2[low, ])),
    theta_high = sqrt(sum(err2[high, ]) / sum(theta2[high, ])),
    # raw and identified beta against the identified truth of these rows,
    # and against the generated beta
    beta_raw_t = rel_err(raw$beta, beta_t),
    beta_id_t = rel_err(id$beta, beta_t),
    beta_raw = rel_err(raw$beta, d$beta),
    beta_id = rel_err(id$beta, d$beta),
    # the same for Gamma, on the test rows
    gamma_raw_t = rel_err(raw$gamma[test, ], gamma_t[test, ]),
    gamma_id_t = rel_err(id$gamma[test, ], gamma_t[test, ]),
    gamma_raw = rel_err(raw$gamma[test, ], d$gamma[test, ]),
    gamma_id = rel_err(id$gamma[test, ], d$gamma[test, ]),
    # raw Gamma on the least and on the most observed rows, against the
    # identified truth
    gamma_low = rel_err(raw$gamma[low, ], gamma_t[low, ]),
    gamma_high = rel_err(raw$gamma[high, ], gamma_t[high, ]),
    # raw against identified estimates
    beta_raw_id = rel_err(raw$beta, id$beta),
    gamma_raw_id = rel_err(raw$gamma[test, ], id$gamma[test, ]),
    lambda_m = out$fit$meta$lambdas[["M"]],
    rank_fit = sum(out$fit$coefficients$d > 0),
    n_iter = out$fit$meta$n_iter, converged = out$fit$meta$converged,
    secs = secs
  ))
}

# One replication ----
# Every size under every mask, on one data set.
run_replication <- function(seed, cfg) {
  d <- make_data(seed, cfg)
  cells <- expand.grid(n = cfg$n, k = seq_along(cfg$kappa_x))
  rows <- lapply(seq_len(nrow(cells)), function(i) {
    return(fit_cell(d, cells$n[i], cells$k[i], seed, cfg))
  })
  return(list(rows = dplyr::bind_rows(rows)))
}

# Run ----
# One task per replication; the seed of a replication is its number.
tasks <- as.list(seq_len(cfg$reps))
names(tasks) <- sprintf("rep_%04d", seq_len(cfg$reps))
task_dir <- file.path(cfg$out_dir, "missingness")
run_tasks(tasks, run_replication, task_dir, cfg$n_workers, cfg = cfg)

results <- list(
  rows = collect_tasks(task_dir)$rows, config = cfg, date = Sys.time(),
  imr_version = as.character(packageVersion("IMR")),
  r_version = R.version.string
)
saveRDS(results, file.path(cfg$out_dir, "results_missingness.rds"))
