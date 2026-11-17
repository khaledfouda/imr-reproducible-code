# Sensitivity analysis, identified components: run the simulations.
#
# Question: can the coefficients that the algorithm returns (raw) be read as
# they are, or only after the identification step (identified)?
#
# The signal is drawn once: X, Z, beta, Gamma, A and B. It has two versions
# with the same X, Z and beta:
#   conditions hold    : the components satisfy (C1)-(C3);
#   conditions violated: the components as generated, which do not.
# T(theta) is the decomposition of theta that satisfies the conditions: the
# generating one in the first version, another one in the second.
#
# One replication draws the noise and the mask, and fits both versions in
# every cell of two scenarios:
#   noise      : the missing rate is fixed and the SNR varies;
#   missingness: the SNR is fixed and the missing rate varies.
# The cells of a replication share the standard normal noise, which the SNR
# scales, and the uniforms behind the mask, so their masks are nested.
#
# Run from the root of the simulation repository, in a terminal, with one
# BLAS thread per worker:
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/identification_run.R lite
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/identification_run.R draft
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/identification_run.R full
# A run that is stopped and started again continues where it stopped: each
# replication has its own result file. Tables and figures:
# identification_results.R.

size <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(size)) size <- "lite"
stopifnot(size %in% c("lite", "draft", "full"))

# Setup ----
# The development version of IMR: the package root is two folders up.
# Replace this line by library(IMR) once that version is installed.
pkgload::load_all("../..", quiet = TRUE)

# generate_simulated_data() is defined in code/simulation/helper.R, which
# also attaches the plotting packages. Only the function is evaluated here.
for (e in parse("./code/simulation/helper.R")) {
  is_generator <- is.call(e) && identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name("generate_simulated_data"))
  if (is_generator) eval(e)
}

# run_tasks() and collect_tasks(): parallel workers, one result file per task
source("./code/sensitivity/runner.R")

# Settings ----
cfg <- list(
  n = 600, m = 600, p = 3, q = 3, r = 4,
  signal_seed = 1,                # seed of X, Z, beta, Gamma, A and B
  snr = c(0.2, 0.5, 1, 3, 10),    # scenario "noise": sd(theta) / sd(noise)
  missing = c(0.1, 0.3, 0.5, 0.7, 0.9),   # scenario "missingness"
  snr_base = 3,                   # SNR of the scenario "missingness"
  missing_base = 0.3,             # missing rate of the scenario "noise"
  # Lasso penalties: 0 is no penalty; c(0, NA, 20) tunes one on 20 values
  beta_grid = 0, gamma_grid = 0,
  val_prop = 0.2, thresh = 1e-6, maxit = 2000,
  n_followed = 6,   # columns of beta, and rows of Gamma, followed one by one
  reps = 100, n_workers = 8,
  out_dir = "./output/sensitivity"
)
if (size == "draft") cfg$reps <- 20
if (size == "lite") {
  cfg$n <- 200
  cfg$m <- 200
  cfg$snr <- c(0.5, 3)
  cfg$missing <- c(0.3, 0.7)
  cfg$reps <- 2
  cfg$n_workers <- 2
  cfg$out_dir <- "./output/sensitivity_lite"
}
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  warning("OPENBLAS_NUM_THREADS is not 1: each worker may use every core.")
}

# The coefficients followed one by one: all of beta in a few columns and all
# of Gamma in a few rows, drawn once. coef_index names them in the order in
# which followed() returns them.
set.seed(cfg$signal_seed)
cfg$cols <- sort(sample(cfg$m, cfg$n_followed))
cfg$rows <- sort(sample(cfg$n, cfg$n_followed))
coef_index <- rbind(
  data.frame(component = "beta",
             covariate = rep(seq_len(cfg$p), times = cfg$n_followed),
             unit = rep(cfg$cols, each = cfg$p)),
  data.frame(component = "gamma",
             covariate = rep(seq_len(cfg$q), each = cfg$n_followed),
             unit = rep(cfg$rows, times = cfg$q))
)
coef_index$coef <- seq_len(nrow(coef_index))

followed <- function(parts, cfg) {
  return(c(parts$beta[, cfg$cols], parts$gamma[cfg$rows, ]))
}

# Relative error ||a - b||_F / ||b||_F
rel_err <- function(a, b) {
  return(sqrt(sum((a - b)^2) / sum(b^2)))
}

# Data ----
# The data set of one cell. The signal depends on cfg$signal_seed and on
# `identified` only; the noise and the mask depend on the replication, and
# have seeds of their own. The generator returns theta, its generating
# components and T(theta) (`identified_truth`).
make_data <- function(rep, snr, missing, identified, cfg) {
  return(generate_simulated_data(
    n = cfg$n, m = cfg$m, r = cfg$r, p = cfg$p, q = cfg$q,
    sparsity = missing, snr = snr, identified = identified,
    seed = cfg$signal_seed, seed_noise = rep, seed_mask = 1e6 + rep
  ))
}

# Distance of a decomposition from (C1)-(C3): the norms of X'M, MZ and
# X'Gamma, each divided by the norms of its two factors. 0: the condition
# holds.
conditions <- function(parts, d) {
  ratio <- function(product, a, b) {
    if (all(product == 0)) return(0)
    return(sqrt(sum(product^2) / (sum(a^2) * sum(b^2))))
  }
  return(c(
    c1 = ratio(crossprod(d$X, parts$M), d$X, parts$M),
    c2 = ratio(parts$M %*% d$Z, parts$M, d$Z),
    c3 = ratio(crossprod(d$X, parts$gamma), d$X, parts$gamma)
  ))
}


# rep = 1; snr = 3; missing = 0.7; identified = FALSE;

# One fit ----
# One cell of one replication: the tuning run, then the raw and the
# identified estimates against the generating components and against
# T(theta). Returns the row of results and the followed coefficients.
fit_cell <- function(rep, snr, missing, identified, cfg) {
  d <- make_data(rep, snr, missing, identified, cfg)
  truth <- d$identified_truth

  # lambda_m is tuned on a validation set; the rank is fixed
  start <- Sys.time()
  out <- suppressMessages(imr(
    d$Y, d$X, d$Z,
    center = FALSE, scale = TRUE, val_prop = cfg$val_prop,
    grid = imr_tune_grid(beta = cfg$beta_grid, gamma = cfg$gamma_grid,
                         rank = cfg$r),
    maxit = cfg$maxit, thresh = cfg$thresh,
    n_cores = 1, seed = rep, verbose = 0
  ))
  secs <- as.numeric(Sys.time() - start, units = "secs")

  # the same fitted matrix, split in two ways
  raw <- reconstruct(out$fit, out$data, trace = FALSE)
  id <- reconstruct(imr_identify(out$fit, out$data), out$data, trace = FALSE)

  err2 <- (raw$estimates - d$theta)^2
  unobs <- d$mask == 0
  row <- data.frame(
    rep = rep, snr = snr, missing = missing, identified = identified,
    obs_rate = mean(d$mask),
    # fitted matrix: all entries, and the unobserved ones
    theta_all = sqrt(sum(err2) / sum(d$theta^2)),
    theta_unobs = sqrt(sum(err2[unobs]) / sum(d$theta[unobs]^2)),
    # check: the identification step leaves the fitted matrix unchanged
    theta_same = rel_err(id$estimates, raw$estimates)
  )

  # <component>_<estimate>_<target>: the raw or identified (id) estimate
  # against the generating component (gen) or its part of T(theta) (T).
  # <component>_raw_id: the distance between the two estimates.
  for (name in c("beta", "gamma", "M")) {
    row[[paste0(name, "_raw_gen")]] <- rel_err(raw[[name]], d[[name]])
    row[[paste0(name, "_id_gen")]] <- rel_err(id[[name]], d[[name]])
    row[[paste0(name, "_raw_T")]] <- rel_err(raw[[name]], truth[[name]])
    row[[paste0(name, "_id_T")]] <- rel_err(id[[name]], truth[[name]])
    row[[paste0(name, "_raw_id")]] <- rel_err(raw[[name]], id[[name]])
  }

  # c1_raw, ..., c3_id: distance of each split from (C1)-(C3)
  row[c("c1_raw", "c2_raw", "c3_raw")] <- as.list(conditions(raw, d))
  row[c("c1_id", "c2_id", "c3_id")] <- as.list(conditions(id, d))

  # Shares of the squared error of the fitted matrix that the identified
  # beta, Gamma and M carry. The three parts are orthogonal: they sum to 1.
  row$share_beta <- sum((d$X %*% (id$beta - truth$beta))^2) / sum(err2)
  row$share_gamma <- sum(tcrossprod(id$gamma - truth$gamma, d$Z)^2) /
    sum(err2)
  row$share_M <- sum((id$M - truth$M)^2) / sum(err2)

  row$lambda_m <- out$fit$meta$lambdas[["M"]]
  row$rank_fit <- sum(out$fit$coefficients$d > 0)
  row$n_iter <- out$fit$meta$n_iter
  row$converged <- out$fit$meta$converged
  row$secs <- secs

  coefs <- data.frame(
    rep = rep, snr = snr, missing = missing, identified = identified,
    coef = coef_index$coef, raw = followed(raw, cfg), id = followed(id, cfg)
  )
  return(list(row = row, coefs = coefs))
}

# One replication ----
# Both versions of the truth in every cell of the two scenarios. The cell
# (snr_base, missing_base) belongs to both scenarios and is fitted once.
run_replication <- function(rep, cfg) {
  cells <- unique(rbind(
    data.frame(snr = cfg$snr, missing = cfg$missing_base),
    data.frame(snr = cfg$snr_base, missing = cfg$missing)
  ))
  # no common column: every cell with each version of the truth
  cells <- merge(cells, data.frame(identified = c(TRUE, FALSE)))
  fits <- lapply(seq_len(nrow(cells)), function(i) {
    return(fit_cell(rep, cells$snr[i], cells$missing[i],
                    cells$identified[i], cfg))
  })
  return(list(
    rows = dplyr::bind_rows(lapply(fits, `[[`, "row")),
    extra = dplyr::bind_rows(lapply(fits, `[[`, "coefs"))
  ))
}

# Truth ----
# What no fit is needed for, in each version of the truth: the gap between
# every generating component and its part of T(theta), and the two values of
# the followed coefficients.
describe_truth <- function(identified, cfg) {
  d <- make_data(1, cfg$snr_base, cfg$missing_base, identified, cfg)
  truth <- d$identified_truth
  return(list(
    facts = data.frame(
      identified = identified, sd_theta = sd(d$theta), rank_theta = d$rank,
      gap_beta = rel_err(truth$beta, d$beta),
      gap_gamma = rel_err(truth$gamma, d$gamma),
      gap_M = rel_err(truth$M, d$M)
    ),
    coefs = cbind(
      identified = identified, coef_index,
      gen = followed(d, cfg), T_theta = followed(truth, cfg)
    )
  ))
}

# Run ----
# One task per replication; the seed of a replication is its number.
tasks <- as.list(seq_len(cfg$reps))
names(tasks) <- sprintf("rep_%04d", seq_len(cfg$reps))
task_dir <- file.path(cfg$out_dir, "identification")
run_tasks(tasks, run_replication, task_dir, cfg$n_workers, cfg = cfg)

collected <- collect_tasks(task_dir)
truths <- lapply(c(TRUE, FALSE), describe_truth, cfg = cfg)
results <- list(
  rows = collected$rows,
  coefs = dplyr::bind_rows(collected$extra),
  truth = list(
    facts = dplyr::bind_rows(lapply(truths, `[[`, "facts")),
    coefs = dplyr::bind_rows(lapply(truths, `[[`, "coefs"))
  ),
  config = cfg, date = Sys.time(),
  imr_version = as.character(packageVersion("IMR")),
  r_version = R.version.string
)
saveRDS(results, file.path(cfg$out_dir, "results_identification.rds"))
