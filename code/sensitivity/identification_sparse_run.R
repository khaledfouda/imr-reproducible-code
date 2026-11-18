# Sensitivity analysis, identified components with sparse coefficients: run
# the simulations.
#
# The simulations of identification_run.R, with two changes:
#   - half of the entries of beta and half of the entries of Gamma are zero,
#     at positions drawn at random, once for each signal;
#   - lambda_beta and lambda_gamma are tuned on the validation set, together
#     with lambda_m.
# Only the truth as generated is used. The projection that makes a truth
# satisfy (C3) would remove the zeros of Gamma.
#
# Everything that is not defined below comes from identification_run.R: the
# setup, the settings (dimensions, signals, replications, scenarios, sizes
# of the run), the followed coefficients, rel_err(), conditions() and
# describe_truth().
#
# Run from the root of the simulation repository, in a terminal, with one
# BLAS thread per worker and the size lite, draft or full:
#   OPENBLAS_NUM_THREADS=1 Rscript \
#     code/sensitivity/identification_sparse_run.R lite
# A run that is stopped and started again continues where it stopped: each
# replication of each signal has its own result file. Tables and figures:
# identification_sparse_results.R.

# Shared code ----
# The setup, the settings and the functions of identification_run.R: every
# expression of that file before its section "Run", which starts by
# defining `tasks`. It reads the size of the run from the command line.
shared <- parse("./code/sensitivity/identification_run.R")
is_run <- vapply(shared, function(e) {
  return(is.call(e) && identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name("tasks")))
}, logical(1))
for (e in shared[seq_len(which(is_run)[1] - 1)]) eval(e)

# Settings ----
# Added to, or changed in, the settings of identification_run.R.
cfg$zeros <- 0.5   # share of the entries of beta, and of Gamma, that are zero
# Grids of lambda_beta and lambda_gamma: lasso_length values from 0 to
# lasso_top times the noise level of the gradient (see fit_cell()). The
# automatic grid of the package goes up to the value that sets every
# coefficient to zero, which is too coarse near the useful values.
cfg$lasso_top <- 2
cfg$lasso_length <- 10
if (size == "lite") cfg$lasso_length <- 5

# Data ----
# The data set of one cell, with sparse coefficients. The generator gives
# X, Z, beta, Gamma, A, B and the mask of identification_run.R. The zeros
# are then set, and theta, the noise and T(theta) are computed from the
# sparse coefficients as the generator computes them. `identified` is kept
# as an argument because describe_truth() passes it.
make_data <- function(signal, draw, snr, missing, identified, cfg) {
  stopifnot(!identified)
  d <- generate_simulated_data(
    n = cfg$n, m = cfg$m, r = cfg$r, p = cfg$p, q = cfg$q,
    sparsity = missing, snr = Inf, identified = FALSE,
    seed = signal, seed_mask = 1e6 + draw, centered = TRUE
  )

  # the zeros: positions drawn once for the signal
  zero <- withr::with_seed(2e6 + signal, list(
    beta = sample(length(d$beta), round(cfg$zeros * length(d$beta))),
    gamma = sample(length(d$gamma), round(cfg$zeros * length(d$gamma)))
  ))
  d$beta[zero$beta] <- 0
  d$gamma[zero$gamma] <- 0
  d$theta <- d$X %*% d$beta + tcrossprod(d$gamma, d$Z) + d$M

  # noise: the standard normal draws of the generator for this seed, with
  # standard deviation sd(theta) / snr
  noise <- withr::with_seed(draw, matrix(rnorm(cfg$n * cfg$m), cfg$n, cfg$m))
  d$Y <- (d$theta + (sd(d$theta) / snr) * noise) * d$mask
  d$params$snr <- snr

  # T(theta): beta = (X'X)^{-1} X' theta, Gamma = Px_perp theta Z (Z'Z)^{-1}
  # and M = Px_perp theta Pz_perp. It is not sparse.
  qx <- qr(d$X)
  qz <- qr(d$Z)
  rest <- qr.resid(qx, d$theta)
  d$identified_truth <- list(
    M = t(qr.resid(qz, t(rest))), beta = qr.coef(qx, d$theta),
    gamma = t(qr.coef(qz, t(rest)))
  )
  return(d)
}

# One fit ----
# fit_cell() of identification_run.R, with lambda_beta and lambda_gamma
# tuned, and with what the Lasso adds to the results: the selected
# penalties and the zeros of the raw estimates.
fit_cell <- function(signal, rep, snr, missing, identified, cfg) {
  # one number per signal and replication: the seed of the noise, of the
  # mask and of the validation split
  draw <- 1000 * (signal - 1) + rep
  d <- make_data(signal, draw, snr, missing, identified, cfg)
  truth <- d$identified_truth

  # The tuned fit of identification_run.R: lambda_m tuned on a validation
  # set, the rank fixed, the covariates already centered in the data.
  # `lasso` is c(min, max, length) for lambda_beta and for lambda_gamma.
  tuned_fit <- function(lasso_beta, lasso_gamma) {
    return(suppressMessages(imr(
      d$Y, d$X, d$Z,
      center = FALSE, scale = TRUE, val_prop = cfg$val_prop,
      grid = imr_tune_grid(beta = lasso_beta, gamma = lasso_gamma,
        rank = cfg$r),
      maxit = cfg$maxit, thresh = cfg$thresh,
      n_cores = 1, seed = draw, verbose = 0
    )))
  }
  start <- Sys.time()

  # First without Lasso penalties. Its residuals give the noise level, and
  # its raw estimates are kept as the reference "no Lasso".
  plain <- tuned_fit(0, 0)
  raw_plain <- reconstruct(plain$fit, plain$data, trace = FALSE)
  sigma_hat <- sqrt(mean(plain$fit$residuals@x^2))

  # Then with lambda_beta and lambda_gamma tuned as well. For a coefficient
  # that is zero, the gradient of the loss is noise with standard deviation
  # sigma * sqrt(number of observed entries of its column, or row), on the
  # scaled covariates. The grids go from 0 to lasso_top times that.
  top_beta <- cfg$lasso_top * sigma_hat * sqrt(mean(colSums(d$mask)))
  top_gamma <- cfg$lasso_top * sigma_hat * sqrt(mean(rowSums(d$mask)))
  out <- tuned_fit(c(0, top_beta, cfg$lasso_length),
    c(0, top_gamma, cfg$lasso_length))
  secs <- as.numeric(Sys.time() - start, units = "secs")

  # the same fitted matrix, split in two ways
  raw <- reconstruct(out$fit, out$data, trace = FALSE)
  id <- reconstruct(imr_identify(out$fit, out$data), out$data, trace = FALSE)

  unobs <- d$mask == 0
  # squared error of the fitted matrix
  total <- sum((raw$estimates - d$theta)^2)
  row <- data.frame(
    signal = signal, rep = rep, snr = snr, missing = missing,
    identified = identified, obs_rate = mean(d$mask),
    # fitted matrix: all entries, and the unobserved ones
    theta_all = rel_err(raw$estimates, d$theta),
    theta_unobs = rel_err(raw$estimates[unobs], d$theta[unobs]),
    # 1 - R^2 of the fitted matrix: its squared error over the variation of
    # theta around its mean
    one_minus_r2 = total / sum((d$theta - mean(d$theta))^2),
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
  row$share_beta <- sum((d$X %*% (id$beta - truth$beta))^2) / total
  row$share_gamma <- sum(tcrossprod(id$gamma - truth$gamma, d$Z)^2) / total
  row$share_M <- sum((id$M - truth$M)^2) / total

  # The zeros of the raw estimates of beta and Gamma:
  #   <name>_zeros        : share of their entries that are exactly zero;
  #   <name>_zeros_found  : share of the true zeros estimated as zero;
  #   <name>_nonzeros_kept: share of the true non-zeros estimated as
  #                         non-zero.
  for (name in c("beta", "gamma")) {
    zero_true <- d[[name]] == 0
    zero_raw <- raw[[name]] == 0
    row[[paste0(name, "_zeros")]] <- mean(zero_raw)
    row[[paste0(name, "_zeros_found")]] <- mean(zero_raw[zero_true])
    row[[paste0(name, "_nonzeros_kept")]] <- mean(!zero_raw[!zero_true])
  }

  # the fit without Lasso penalties, on the same data
  row$theta_all_plain <- rel_err(raw_plain$estimates, d$theta)
  row$beta_plain_gen <- rel_err(raw_plain$beta, d$beta)
  row$gamma_plain_gen <- rel_err(raw_plain$gamma, d$gamma)

  # the selected penalties, the top of the grid of each Lasso penalty and
  # the estimated noise level
  row$lambda_beta <- out$fit$meta$lambdas[["beta"]]
  row$lambda_gamma <- out$fit$meta$lambdas[["gamma"]]
  row$lambda_beta_max <- top_beta
  row$lambda_gamma_max <- top_gamma
  row$sigma_hat <- sigma_hat
  row$sigma <- sd(d$theta) / snr
  row$lambda_m <- out$fit$meta$lambdas[["M"]]
  row$rank_fit <- sum(out$fit$coefficients$d > 0)
  row$n_iter <- out$fit$meta$n_iter
  row$converged <- out$fit$meta$converged
  row$secs <- secs

  coefs <- data.frame(
    signal = signal, rep = rep, snr = snr, missing = missing,
    identified = identified, coef = coef_index$coef,
    raw = followed(raw, cfg), id = followed(id, cfg)
  )
  return(list(row = row, coefs = coefs))
}

# One replication ----
# `task` names a signal and one of its replications. The truth as generated
# is fitted in every cell of the two scenarios. The cell (snr_base,
# missing_base) belongs to both scenarios and is fitted once.
run_replication <- function(task, cfg) {
  cells <- unique(rbind(
    data.frame(snr = cfg$snr, missing = cfg$missing_base),
    data.frame(snr = cfg$snr_base, missing = cfg$missing)
  ))
  fits <- lapply(seq_len(nrow(cells)), function(i) {
    return(fit_cell(task$signal, task$rep, cells$snr[i], cells$missing[i],
      FALSE, cfg))
  })
  return(list(
    rows = dplyr::bind_rows(lapply(fits, `[[`, "row")),
    extra = dplyr::bind_rows(lapply(fits, `[[`, "coefs"))
  ))
}

# Run ----
# One task per signal and replication.
tasks <- expand.grid(rep = seq_len(cfg$reps), signal = cfg$signals)
tasks <- split(tasks, sprintf("signal_%02d_rep_%02d", tasks$signal, tasks$rep))
task_dir <- file.path(cfg$out_dir, "identification_sparse")
run_tasks(tasks, run_replication, task_dir, cfg$n_workers, cfg = cfg)

# The results of these tasks only: the folder may hold the files of a run
# with other settings.
files <- file.path(task_dir, paste0(names(tasks), ".rds"))
done <- lapply(files[file.exists(files)], readRDS)
truths <- lapply(cfg$signals, describe_truth, identified = FALSE, cfg = cfg)
results <- list(
  rows = dplyr::bind_rows(lapply(done, `[[`, "rows")),
  coefs = dplyr::bind_rows(lapply(done, `[[`, "extra")),
  truth = list(
    facts = dplyr::bind_rows(lapply(truths, `[[`, "facts")),
    coefs = dplyr::bind_rows(lapply(truths, `[[`, "coefs"))
  ),
  config = cfg, date = Sys.time(),
  imr_version = as.character(packageVersion("IMR")),
  r_version = R.version.string
)
saveRDS(results, file.path(cfg$out_dir, "results_identification_sparse.rds"))
