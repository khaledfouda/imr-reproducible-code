library(tidyverse)
# library(IMR)
library(kableExtra)
library(magrittr)
library(scales)
library(RSSthemes)

#------------------------------------------------------------------------------


source("./code/simulation/config_default.R")
if (file.exists("./code/simulation/config.R")) {
  source("./code/simulation/config.R")
}

generate_simulated_data <- function(
  n = 300,
  m = 400,
  r = 10,
  p = 6,
  q = 6,
  sparsity = 0.8,
  sparsity_beta = 0,
  sparsity_gamma = 0,
  shared = FALSE,
  snr = 1,
  mechanism = c("mcar", "mar"),
  mar_strength = 0.25,
  identified = FALSE,
  row_intercept = FALSE,
  col_intercept = FALSE,
  seed = NULL,
  seed_mask = NULL,
  seed_noise = NULL
) {
  mechanism <- match.arg(mechanism)

  #  checks ---------------------------------------------------
  if (sparsity <= 0 || sparsity >= 1) stop("`sparsity` must be in (0, 1).")
  if (sparsity_beta < 0 || sparsity_beta > 1) stop("`sparsity_beta` must be in [0, 1].")
  if (sparsity_gamma < 0 || sparsity_gamma > 1) stop("`sparsity_gamma` must be in [0, 1].")
  if (snr <= 0) stop("`snr` must be positive. `snr = Inf` gives no noise.")
  if (mar_strength < 0) stop("`mar_strength` must be nonnegative.")
  if (mechanism == "mar" && p == 0 && q == 0) stop("`mechanism = \"mar\"` needs covariates.")
  if (identified && shared) stop("`identified = TRUE` does not support shared effects.")

  # Random numbers -----------------------------------------------------------
  # Each component has its own seed, so two calls that differ in one argument
  # share every other component. The caller's RNG state is left as it was,
  # except for one draw when `seed` is NULL.
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  withr::local_preserve_seed()
  set.seed(seed)
  components <- c(
    "x", "z", "beta", "gamma", "zero_beta", "zero_gamma",
    "beta0", "gamma0", "a", "b", "prob", "mask", "noise"
  )
  seeds <- sample.int(.Machine$integer.max, length(components))
  names(seeds) <- components
  if (!is.null(seed_mask)) seeds["mask"] <- seed_mask
  if (!is.null(seed_noise)) seeds["noise"] <- seed_noise

  # Covariates X and Z -------------------------------------------------------
  set.seed(seeds[["x"]])
  x_mat <- if (p > 0) matrix(runif(n * p), nrow = n, ncol = p) else NULL
  set.seed(seeds[["z"]])
  z_mat <- if (q > 0) matrix(runif(m * q), nrow = m, ncol = q) else NULL

  # Coefficient matrices -----------------------------------------------------
  # Row k of beta is N(mean_k, sd_k^2), with |mean_k| ~ U(0.1, 1) and
  # sd_k ~ U(0.5, 1). Column k of gamma likewise. A shared effect is the
  # vector of means.
  beta_mat <- NULL
  if (p > 0) {
    set.seed(seeds[["beta"]])
    beta_means <- runif(p, 0.1, 1) * sample(c(-1, 1), p, replace = TRUE)
    if (shared) {
      beta_mat <- matrix(beta_means, p, 1)
    } else {
      beta_sds <- runif(p, 0.5, 1)
      beta_mat <- beta_means + beta_sds * matrix(rnorm(p * m), nrow = p, ncol = m)
    }
  }

  gamma_mat <- NULL
  if (q > 0) {
    set.seed(seeds[["gamma"]])
    gamma_means <- runif(q, 0.1, 1) * sample(c(-1, 1), q, replace = TRUE)
    if (shared) {
      gamma_mat <- matrix(gamma_means, q, 1)
    } else {
      gamma_sds <- runif(q, 0.5, 1)
      gamma_mat <- t(gamma_means + gamma_sds * matrix(rnorm(q * n), nrow = q, ncol = n))
    }
  }

  # Zeros in the coefficients (set before any projection) --------------------
  if (sparsity_beta > 0 && p > 0) {
    set.seed(seeds[["zero_beta"]])
    to_zero_b <- sample(length(beta_mat), size = round(sparsity_beta * length(beta_mat)))
    beta_mat[to_zero_b] <- 0
  }

  if (sparsity_gamma > 0 && q > 0) {
    set.seed(seeds[["zero_gamma"]])
    to_zero_g <- sample(length(gamma_mat), size = round(sparsity_gamma * length(gamma_mat)))
    gamma_mat[to_zero_g] <- 0
  }

  # Intercepts: beta0 for the rows, gamma0 for the columns -------------------
  beta0 <- gamma0 <- NULL
  if (row_intercept) {
    set.seed(seeds[["beta0"]])
    beta0 <- rnorm(n)
  }
  if (col_intercept) {
    set.seed(seeds[["gamma0"]])
    gamma0 <- rnorm(m)
  }

  # Low-rank factors A and B, with M = A B' ----------------------------------
  set.seed(seeds[["a"]])
  a_mat <- matrix(runif(n * r, -1, 1), nrow = n, ncol = r)
  set.seed(seeds[["b"]])
  b_mat <- matrix(runif(m * r, -1, 1), nrow = m, ncol = r)

  # Conditions (C1)-(C3) -----------------------------------------------------
  # xt = [1_n, X] when the column intercepts are present, and zt = [1_m, Z]
  # when the row intercepts are, as in IMR::imr_identify(). With
  # `identified = TRUE`, A, B and [beta0, gamma] are projected so that
  # xt' A = 0, zt' B = 0 and xt' [beta0, gamma] = 0. The projection of gamma
  # removes its zeros.
  xt <- cbind(if (col_intercept) rep(1, n), x_mat)
  zt <- cbind(if (row_intercept) rep(1, m), z_mat)
  qx <- if (!is.null(xt)) qr(xt)
  qz <- if (!is.null(zt)) qr(zt)
  if (identified) {
    if (!is.null(qx)) {
      a_mat <- qr.resid(qx, a_mat)
      if (q > 0) gamma_mat <- qr.resid(qx, gamma_mat)
      if (row_intercept) beta0 <- qr.resid(qx, beta0)
    }
    if (!is.null(qz)) b_mat <- qr.resid(qz, b_mat)
  }
  m_mat <- tcrossprod(a_mat, b_mat)

  # Theta --------------------------------------------------------------------
  theta <- m_mat
  if (p > 0 && !shared) theta <- theta + (x_mat %*% beta_mat)
  if (q > 0 && !shared) theta <- theta + (gamma_mat %*% t(z_mat))
  if (p > 0 && shared) theta <- sweep(theta, 1, as.vector(x_mat %*% beta_mat), "+")
  if (q > 0 && shared) theta <- sweep(theta, 2, as.vector(gamma_mat %*% t(z_mat)), "+")
  if (row_intercept) theta <- theta + beta0
  if (col_intercept) theta <- sweep(theta, 2, gamma0, "+")

  # Noise --------------------------------------------------------------------
  # Its standard deviation is sd(theta) / snr, so snr = Inf gives no noise.
  y_mat <- theta
  if (is.finite(snr)) {
    set.seed(seeds[["noise"]])
    y_mat <- y_mat + (sd(theta) / snr) * matrix(rnorm(n * m), nrow = n, ncol = m)
  }

  # Missingness mask (1 = observed, 0 = missing) -----------------------------
  # Entry (i, j) is observed when u_ij < prob_ij. The uniforms u are the same
  # for every mechanism and rate, so the masks are nested in `sparsity`.
  # MCAR: prob = 1 - sparsity.
  # MAR: logit(prob_ij) = a0 + mar_strength * (row_score_i + col_score_j),
  # where the scores are the standardized covariates times a random unit
  # vector, and a0 is set so that mean(prob) = 1 - sparsity.
  prob <- 1 - sparsity
  if (mechanism == "mar") {
    set.seed(seeds[["prob"]])
    row_score <- numeric(n)
    col_score <- numeric(m)
    if (p > 0) {
      w_x <- rnorm(p)
      row_score <- drop(scale(x_mat) %*% (w_x / sqrt(sum(w_x^2))))
    }
    if (q > 0) {
      w_z <- rnorm(q)
      col_score <- drop(scale(z_mat) %*% (w_z / sqrt(sum(w_z^2))))
    }
    score <- mar_strength * outer(row_score, col_score, "+")
    a0 <- uniroot(
      function(a) mean(plogis(a + score)) - (1 - sparsity),
      interval = c(-10, 10), extendInt = "upX", tol = 1e-10
    )$root
    prob <- plogis(a0 + score)
  }
  set.seed(seeds[["mask"]])
  mask <- matrix((runif(n * m) < prob) * 1L, nrow = n, ncol = m)
  if (min(rowSums(mask), colSums(mask)) < 2) {
    stop("A row or column has fewer than two observed entries: lower `sparsity` or `mar_strength`.")
  }
  y_mat <- y_mat * mask

  # Identified truth T(theta) ------------------------------------------------
  # The decomposition of theta that satisfies (C1)-(C3) for xt and zt:
  #   [gamma0'; beta] = (xt' xt)^{-1} xt' theta,
  #   [beta0, gamma]  = Px_perp theta zt (zt' zt)^{-1},
  #   M               = Px_perp theta Pz_perp.
  # It equals the generating components when `identified = TRUE`.
  identified_truth <- NULL
  if (!shared) {
    b_t <- g_t <- NULL
    m_t <- theta
    if (!is.null(qx)) {
      b_t <- qr.coef(qx, theta)
      m_t <- qr.resid(qx, theta)
    }
    if (!is.null(qz)) {
      g_t <- t(qr.coef(qz, t(m_t)))
      m_t <- t(qr.resid(qz, t(m_t)))
    }
    identified_truth <- list(
      M      = m_t,
      beta   = if (p > 0) b_t[col_intercept + seq_len(p), , drop = FALSE],
      gamma  = if (q > 0) g_t[, row_intercept + seq_len(q), drop = FALSE],
      beta0  = if (row_intercept) g_t[, 1],
      gamma0 = if (col_intercept) b_t[1, ]
    )
  }

  # Rank of theta ------------------------------------------------------------
  # The columns of `left` span the column space of theta and those of `right`
  # its row space, so theta and the small matrix Ql' theta Qr have the same rank.
  left <- cbind(1, x_mat, if (!shared) gamma_mat, a_mat, beta0)
  right <- cbind(1, z_mat, if (!shared && p > 0) t(beta_mat), b_mat, gamma0)
  rank_theta <- qr(crossprod(qr.Q(qr(left)), theta) %*% qr.Q(qr(right)))$rank

  #  Done :) ------------------------------------------------------------------
  out <- list(
    theta = theta,
    mask  = mask,
    Y     = y_mat,
    M     = m_mat,
    A     = a_mat,
    B     = b_mat,
    rank  = rank_theta,
    prob  = prob,
    identified_truth = identified_truth,
    params = mget(names(formals()))
  )

  if (p > 0) {
    out$X <- x_mat
    out$beta <- beta_mat
  }
  if (q > 0) {
    out$Z <- z_mat
    out$gamma <- gamma_mat
  }
  out$beta0 <- beta0
  out$gamma0 <- gamma0

  out
}
