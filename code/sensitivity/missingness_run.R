# Sensitivity analysis, missingness mechanism: run the simulations.
#
# One replication:
#   1. a complete data set with n_max rows, theta = X beta + A B' and
#      Y = theta + noise, where X'A = 0 on the n_max rows;
#   2. four masks on it, MCAR and three MAR mechanisms of increasing
#      strength, all drawn from the same uniforms;
#   3. for each mask and each n, a tuned fit on the first n rows.
# The fits of a replication share everything except the mask and n. The
# error of the fitted matrix is always measured on the first min(n) rows.
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
  m = 400, p = 4, r = 4, snr = 3,
  rate = 0.3,                   # mean probability that an entry is observed
  kappa = c(0, 0.5, 0.75, 1),   # strength of the mechanism; 0 is MCAR
  beta_grid = c(0, NA, 20),     # lambda_beta: min, max (NA: automatic), length
  val_prop = 0.2, thresh = 1e-6, maxit = 2000,
  reps = 500, n_workers = 8,
  out_dir = "./output/sensitivity"
)
if (size == "draft") cfg$reps <- 27
if (size == "lite") {
  cfg$n <- c(100, 200, 300)
  cfg$m <- 200
  cfg$beta_grid <- c(0, NA, 4)
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
# The complete data set of a replication: X, beta, theta and Y with n_max
# rows, the uniforms behind every mask, and the probability that an entry of
# each row is observed (one column per kappa).
make_data <- function(seed, cfg) {
  n_max <- max(cfg$n)
  # centered: each column of X has mean zero on the n_max rows.
  # identified: A is orthogonal to X on the n_max rows, condition (C1).
  # On the first n < n_max rows, both hold only approximately.
  d <- generate_simulated_data(
    n = n_max, m = cfg$m, r = cfg$r, p = cfg$p, q = 0, seed = seed,
    centered = TRUE, identified = TRUE
  )
  theta <- d$theta   # X beta + A B'

  # noise, direction of the MAR score and uniforms, with a seed of their own
  set.seed(1e6 + seed)
  sigma <- sd(theta) / cfg$snr
  noise <- matrix(rnorm(n_max * cfg$m, sd = sigma), n_max, cfg$m)
  w <- rnorm(cfg$p)
  u <- matrix(runif(n_max * cfg$m), n_max, cfg$m)

  # MAR: logit(prob_i) = a0 + kappa * score_i, where the score is the
  # standardized X times a random unit vector and a0 gives mean(prob) = rate.
  # kappa = 0 gives prob = rate in every row (MCAR).
  score <- drop(scale(d$X) %*% (w / sqrt(sum(w^2))))
  prob <- vapply(cfg$kappa, function(kappa) {
    a0 <- uniroot(
      function(a) mean(plogis(a + kappa * score)) - cfg$rate,
      interval = c(-10, 10), tol = 1e-10
    )$root
    return(plogis(a0 + kappa * score))
  }, numeric(n_max))

  return(list(
    X = d$X, beta = d$beta, theta = theta, Y = theta + noise,
    u = u, prob = prob
  ))
}

# One fit ----
# The first n rows under the mask of cfg$kappa[k]: the tuning run and the
# errors, as one row of results.
fit_cell <- function(d, n, k, seed, cfg) {
  rows <- seq_len(n)
  test <- seq_len(min(cfg$n))   # the rows that every size contains
  x_mat <- d$X[rows, ]
  # the first n rows of the mask of the complete data set (TRUE = observed),
  # so the masks of the different sizes are nested
  mask <- d$u[rows, ] < d$prob[rows, k]
  cell <- data.frame(
    seed = seed, n = n, kappa = cfg$kappa[k],
    pi_ratio = max(d$prob[rows, k]) / min(d$prob[rows, k]),
    obs_rate = mean(mask),
    min_row_obs = min(rowSums(mask)), min_col_obs = min(colSums(mask)),
    # distance between beta and the identified truth of these rows,
    # (X'X)^{-1} X' theta: zero when X'A = 0 on them, as at n = n_max
    beta_gap = rel_err(qr.coef(qr(x_mat), d$theta[rows, ]), d$beta)
  )
  # imr() needs two observed entries in every row and column
  if (min(cell$min_row_obs, cell$min_col_obs) < 2) {
    return(cbind(cell, status = "skipped"))
  }

  # lambda_beta and lambda_m are tuned on a validation set; the rank is fixed
  start <- Sys.time()
  out <- suppressMessages(imr(
    d$Y[rows, ] * mask, x_mat,
    center = FALSE, scale = TRUE, val_prop = cfg$val_prop,
    grid = imr_tune_grid(beta = cfg$beta_grid, rank = cfg$r),
    maxit = cfg$maxit, thresh = cfg$thresh,
    n_cores = 1, seed = seed, verbose = 0
  ))
  secs <- as.numeric(Sys.time() - start, units = "secs")

  # raw estimates, and identified ones: beta tilde = (X'X)^{-1} X' theta hat
  raw <- reconstruct(out$fit, out$data, trace = FALSE)
  id <- reconstruct(imr_identify(out$fit, out$data), out$data, trace = FALSE)

  err2 <- (raw$estimates[test, ] - d$theta[test, ])^2
  theta2 <- d$theta[test, ]^2
  unobs <- !mask[test, ]
  return(cbind(
    cell, status = "ok",
    # fitted matrix on the test rows: all entries, and the unobserved ones
    theta_all = sqrt(sum(err2) / sum(theta2)),
    theta_unobs = sqrt(sum(err2[unobs]) / sum(theta2[unobs])),
    # raw and identified coefficients against the true beta, and each other
    beta_raw = rel_err(raw$beta, d$beta),
    beta_id = rel_err(id$beta, d$beta),
    beta_raw_id = rel_err(raw$beta, id$beta),
    lambda_beta = out$fit$meta$lambdas[["beta"]],
    lambda_beta_max = out$grid$beta$max,
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
  cells <- expand.grid(n = cfg$n, k = seq_along(cfg$kappa))
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
