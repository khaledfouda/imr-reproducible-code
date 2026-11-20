# Sensitivity analysis, initialization: run the simulations.
#
# Question: with the data and every parameter of the fit fixed, does the fit
# depend on its starting values?
#
# A signal is X, Z, beta, Gamma, A and B, drawn once, with one noise matrix
# and one mask. The covariates are centered and the model has no intercepts.
# A signal has two versions with the same X, Z, beta, noise draws and mask:
#   conditions hold    : the components satisfy (C1)-(C3);
#   conditions violated: the components as generated, which do not.
#
# Nothing is tuned. Every fit is imr_fit() with the true rank, lambda_beta =
# lambda_gamma = 0 and the same lambda_m. In each version of a signal, only
# the starting values change:
#   default: the coefficients of the fit without M, then the truncated SVD of
#            its residuals. It is the start of imr_fit(), built here;
#   zero   : zero coefficients, and the truncated SVD of the zero-filled Y;
#   random : random coefficients and random subspaces, several draws.
# Every start is run at two thresholds: "tight", to compare the limits, and
# "short", the threshold of the other parts.
#
# Along a fit, the objective is saved after every iteration, and the errors
# at a set of iterations (the checkpoints).
#
# Run from the root of the simulation repository, in a terminal, with one
# BLAS thread per worker:
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/initialization_run.R lite
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/initialization_run.R draft
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/initialization_run.R full
# A run that is stopped and started again continues where it stopped: each
# version of each signal has its own result file. Tables and figures:
# initialization_results.R.

size <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(size)) size <- "lite"
stopifnot(size %in% c("lite", "draft", "full"))

# Setup ----
library(IMR)
source("./code/simulation/simulation_function.R")
# run_tasks(): parallel workers, one result file per task
source("./code/sensitivity/runner.R")

# Settings ----
cfg <- list(
  n = 600, m = 600, p = 3, q = 3, r = 4,
  snr = 3, missing = 0.7,
  # Fixed, not tuned. 0.5 is in the middle of the values that the tuning
  # selects in this design with the rank fixed at r (identification_run.R at
  # 70% missing: median 0.48).
  lambda_m = 0.5,
  signals = 1:10,                 # seeds of the signals
  n_random = 25,                  # random starts of each version of a signal
  thresh = c(tight = 1e-12, short = 1e-6),
  maxit = c(tight = 10000, short = 2000),
  n_checkpoints = 80,   # iterations with saved errors: about this many
  increase_tol = 1e-12,   # a relative increase of the objective above this
  n_workers = 8,
  out_dir = "./output/sensitivity"
)
# draft: the first signals and the first random starts of the full run
if (size == "draft") {
  cfg$signals <- 1:5
  cfg$n_random <- 10
}
if (size == "lite") {
  cfg$n <- 200
  cfg$m <- 200
  cfg$signals <- 1:2
  cfg$n_random <- 2
  cfg$n_workers <- 2
  cfg$out_dir <- "./output/sensitivity_lite"
}
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  warning("OPENBLAS_NUM_THREADS is not 1: each worker may use every core.")
}

# The checkpoints: iterations spread evenly on the log scale, from the first
# to the largest number of iterations.
cfg$checkpoints <- unique(round(
  10^seq(0, log10(max(cfg$maxit)), length.out = cfg$n_checkpoints)
))

# The starts of each version of a signal. The default start comes first: the
# other fits are compared with its fit.
starts <- data.frame(
  kind = c("default", "zero", rep("random", cfg$n_random)),
  start_id = c(0, 0, seq_len(cfg$n_random))
)

# RRMSE, as in the paper: ||a - b||_F / ||b||_F
rel_err <- function(a, b) {
  return(sqrt(sum((a - b)^2) / sum(b^2)))
}

# Data ----
# The data set of one version of a signal: centered covariates, no
# intercepts. The two versions differ in `identified` only, so they share
# X, Z, beta, the standard normal noise and the mask. The generator returns
# theta, its generating components and T(theta) (`identified_truth`).
make_data <- function(signal, identified, cfg) {
  return(generate_simulated_data(
    n = cfg$n, m = cfg$m, r = cfg$r, p = cfg$p, q = cfg$q,
    sparsity = cfg$missing, snr = cfg$snr, identified = identified,
    seed = signal, centered = TRUE
  ))
}

# Starting values ----
# One start, as the list that imr_fit() takes as warm_start: u, d and v with
# r columns, beta and Gamma on the scale of the stored designs, and no
# intercepts. `data` is the imr_data object of the fit.
make_start <- function(kind, start_id, signal, data, level, cfg) {
  n <- nrow(data$Y)
  m <- ncol(data$Y)
  r <- cfg$r
  # one seed per start, the same in the two versions of a signal and at the
  # two thresholds; svd_opt() may draw random numbers
  set.seed(10000 * signal + start_id)
  if (kind == "default") {
    # the fit without M, with the threshold of the fit that follows
    fit0 <- suppressMessages(imr_fit(
      update(data, low_rank_component = FALSE),
      convergence = imr_convergence(
        maxit = cfg$maxit[[level]], thresh = cfg$thresh[[level]]
      )
    ))
    factors <- svd_opt(fit0$residuals, r)
    beta <- fit0$coefficients$beta
    gamma <- fit0$coefficients$gamma
  }
  if (kind == "zero") {
    factors <- svd_opt(data$Y, r)
    beta <- matrix(0, cfg$p, m)
    gamma <- matrix(0, n, cfg$q)
  }
  if (kind == "random") {
    # random subspaces with the leading singular values of the zero-filled
    # Y, and N(0, 1) coefficients
    factors <- list(
      d = svd_opt(data$Y, r)$d,
      u = qr.Q(qr(matrix(rnorm(n * r), n, r))),
      v = qr.Q(qr(matrix(rnorm(m * r), m, r)))
    )
    beta <- matrix(rnorm(cfg$p * m), cfg$p, m)
    gamma <- matrix(rnorm(n * cfg$q), n, cfg$q)
  }
  return(list(
    u = factors$u, d = factors$d, v = factors$v, beta = beta, gamma = gamma,
    beta0 = NULL, gamma0 = NULL
  ))
}

# Errors ----
# The estimates of a fit, in two splits of the same fitted matrix: `raw`, as
# the algorithm returns them, and `id`, after the identification step.
fit_parts <- function(fit, data) {
  return(list(
    raw = reconstruct(fit, data, trace = FALSE),
    id = reconstruct(imr_identify(fit, data), data, trace = FALSE)
  ))
}

# Errors of a fit, as one row: the fitted matrix on all entries and on the
# unobserved ones, the raw beta and Gamma against the generating ones, and
# the identified beta and Gamma against their part of T(theta).
fit_errors <- function(parts, d) {
  unobs <- d$mask == 0
  truth <- d$identified_truth
  return(data.frame(
    theta_all = rel_err(parts$raw$estimates, d$theta),
    theta_unobs = rel_err(parts$raw$estimates[unobs], d$theta[unobs]),
    beta_raw = rel_err(parts$raw$beta, d$beta),
    gamma_raw = rel_err(parts$raw$gamma, d$gamma),
    beta_id = rel_err(parts$id$beta, truth$beta),
    gamma_id = rel_err(parts$id$gamma, truth$gamma)
  ))
}

# Relative differences between the estimates of a fit and those of the
# default-start fit, as one row.
fit_differences <- function(parts, reference) {
  return(data.frame(
    d_theta = rel_err(parts$raw$estimates, reference$raw$estimates),
    d_beta_raw = rel_err(parts$raw$beta, reference$raw$beta),
    d_gamma_raw = rel_err(parts$raw$gamma, reference$raw$gamma),
    d_beta_id = rel_err(parts$id$beta, reference$id$beta),
    d_gamma_id = rel_err(parts$id$gamma, reference$id$gamma)
  ))
}

# One fit ----
# The fit from `start` at one threshold. imr_fit() returns the objective
# after each iteration but not the estimates, so the fit is run in segments
# that end at the checkpoints: each segment starts from the end of the
# previous one, and the errors are computed between two segments. The
# iterations are those of a single call, because imr_fit() keeps nothing
# from one iteration to the next besides the estimates.
# Returns the objective after each iteration, the errors at the checkpoints
# (`path`), the estimates of the last iteration (`parts`) and the facts of
# the fit.
fit_path <- function(data, d, start, level, cfg) {
  maxit <- cfg$maxit[[level]]
  checkpoints <- cfg$checkpoints[cfg$checkpoints <= maxit]
  objective <- numeric(0)
  path <- list()
  done <- 0
  warm_start <- start
  converged <- FALSE
  begin <- Sys.time()
  for (checkpoint in checkpoints) {
    fit <- suppressMessages(imr_fit(
      data, rank = cfg$r, lambda_m = cfg$lambda_m,
      lambda_beta = 0, lambda_gamma = 0,
      convergence = imr_convergence(
        maxit = checkpoint - done, thresh = cfg$thresh[[level]]
      ),
      warm_start = warm_start, return_objective = TRUE
    ))
    # the fit stopped before the checkpoint: it has converged
    converged <- fit$meta$n_iter < checkpoint - done
    done <- done + fit$meta$n_iter
    objective <- c(objective, fit$objective)
    parts <- fit_parts(fit, data)
    path[[length(path) + 1]] <- cbind(iteration = done, fit_errors(parts, d))
    if (converged) break
    warm_start <- fit
  }
  secs <- as.numeric(Sys.time() - begin, units = "secs")

  # relative change of the objective at each iteration: positive = increase
  steps <- diff(objective) / objective[-length(objective)]
  facts <- data.frame(
    n_iter = done, converged = converged,
    objective = objective[length(objective)],
    max_rel_increase = max(steps),
    n_increase = sum(steps > cfg$increase_tol),
    rank_fit = sum(fit$coefficients$d > 0), secs = secs
  )
  return(list(
    objective = objective, path = dplyr::bind_rows(path), parts = parts,
    facts = facts
  ))
}

# One version of a signal ----
# `task` names a signal and a version of its truth. Every start is fitted at
# the two thresholds. Returns one row per fit, the errors at the checkpoints
# of every fit (`history`) and the objectives (`extra`), named after `fit`.
run_version <- function(task, cfg) {
  d <- make_data(task$signal, task$identified, cfg)
  # the covariates are already centered in the data, so the fit does not
  # center them; all the observed entries are used
  data <- update(
    imr_data(d$Y, d$X, d$Z, val_prop = 0, center = FALSE, scale = TRUE),
    row_intercept = FALSE, col_intercept = FALSE
  )

  rows <- list()
  paths <- list()
  objectives <- list()
  for (level in names(cfg$thresh)) {
    reference <- NULL
    for (i in seq_len(nrow(starts))) {
      start <- make_start(starts$kind[i], starts$start_id[i], task$signal,
        data, level, cfg)
      result <- fit_path(data, d, start, level, cfg)
      if (is.null(reference)) reference <- result$parts
      cell <- data.frame(
        fit = sprintf("%s_%s_%02d_%s", task$name, starts$kind[i],
          starts$start_id[i], level),
        signal = task$signal, identified = task$identified,
        kind = starts$kind[i], start_id = starts$start_id[i], level = level
      )
      rows[[cell$fit]] <- cbind(
        cell, thresh = cfg$thresh[[level]], maxit = cfg$maxit[[level]],
        result$facts, fit_errors(result$parts, d),
        fit_differences(result$parts, reference)
      )
      paths[[cell$fit]] <- cbind(cell, result$path)
      objectives[[cell$fit]] <- result$objective
    }
  }
  return(list(
    rows = dplyr::bind_rows(rows), history = dplyr::bind_rows(paths),
    extra = objectives
  ))
}

# Run ----
# One task per signal and version of its truth.
tasks <- expand.grid(identified = c(TRUE, FALSE), signal = cfg$signals)
tasks$name <- sprintf("signal_%02d_%s", tasks$signal,
  ifelse(tasks$identified, "holds", "violated"))
tasks <- split(tasks, tasks$name)
task_dir <- file.path(cfg$out_dir, "initialization")
run_tasks(tasks, run_version, task_dir, cfg$n_workers, cfg = cfg)

# The results of these tasks only: the folder may hold the files of a run
# with other settings.
files <- file.path(task_dir, paste0(names(tasks), ".rds"))
done <- lapply(files[file.exists(files)], readRDS)
results <- list(
  rows = dplyr::bind_rows(lapply(done, `[[`, "rows")),
  paths = dplyr::bind_rows(lapply(done, `[[`, "history")),
  # the objective after each iteration of every fit, named after `fit`
  objective = do.call(c, lapply(done, `[[`, "extra")),
  config = cfg, date = Sys.time(),
  imr_version = as.character(packageVersion("IMR")),
  r_version = R.version.string
)
saveRDS(results, file.path(cfg$out_dir, "results_initialization.rds"))
