# Building blocks shared by all parts: the data set of a cell, the tuned fit,
# and the metrics.
#
# A "cell" is a one-row data frame: the seeds and factor levels of one fit.

# Relative error ||a - b||_F / ||b||_F
rel_err <- function(a, b) {
  return(sqrt(sum((a - b)^2) / sum(b^2)))
}

# Data set of a cell: the baseline, with the arguments of the generator that
# the cell names (n, mechanism, mar_strength, identified, seed, seed_noise,
# seed_mask). A draw with too few observed entries in a row or column is not
# an error here: it returns list(skipped = <message>).
sens_data <- function(cfg, cell) {
  args <- cfg$base
  generator_args <- names(formals(generate_simulated_data))
  for (name in intersect(names(cell), generator_args)) {
    if (!is.na(cell[[name]])) args[[name]] <- cell[[name]]
  }
  d <- tryCatch(
    do.call(generate_simulated_data, args),
    error = function(e) {
      if (!grepl("fewer than two observed", conditionMessage(e))) stop(e)
      return(list(skipped = conditionMessage(e)))
    }
  )
  return(d)
}

# Properties of the data set that need no fit. The gap is the distance from
# the generating component to its part of T(theta).
sens_facts <- function(d) {
  truth <- d$identified_truth
  return(data.frame(
    n_obs = sum(d$mask),
    min_row_obs = min(rowSums(d$mask)),
    min_col_obs = min(colSums(d$mask)),
    pi_min = min(d$prob),
    pi_max = max(d$prob),
    sigma = sd(d$theta) / d$params$snr,
    rank_theta = d$rank,
    gap_beta = rel_err(truth$beta, d$beta),
    gap_gamma = rel_err(truth$gamma, d$gamma),
    gap_M = rel_err(truth$M, d$M)
  ))
}

# One tuning run: lambda_m and the rank by the fast nuclear path, then the
# final fit on all observed entries. `seed_fit` fixes the validation split and
# every random number the fits draw.
sens_tune <- function(d, cfg, seed_fit) {
  set.seed(seed_fit)
  start <- Sys.time()
  out <- suppressMessages(imr(
    d$Y, d$X, d$Z,
    center = cfg$fit$center, scale = cfg$fit$scale,
    val_prop = cfg$fit$val_prop,
    grid = imr_tune_grid(beta = 0, gamma = 0),
    maxit = cfg$fit$maxit, thresh = cfg$fit$thresh,
    n_cores = 1, seed = seed_fit, verbose = 0
  ))
  out$secs <- as.numeric(Sys.time() - start, units = "secs")
  return(out)
}

# Tuning outputs of sens_tune(), as one row
sens_tuning <- function(out) {
  return(data.frame(
    lambda_m = out$params$lambda_m,
    rank_in = out$params$rank_in,       # rank bound of the selected fit
    rank_out = out$params$rank_out,     # its rank on the training entries
    rank_fit = sum(out$fit$coefficients$d > 0),   # rank of the final fit
    valid_error = out$params$verror,
    n_path = nrow(out$history),         # number of fits along the path
    lambda_max = out$grid$nuclear$max,
    n_iter = out$fit$meta$n_iter,
    converged = out$fit$meta$converged,
    time_total = out$secs,
    time_final = out$fit$meta$time_secs
  ))
}

# Estimates of a fit, raw and identified: each is a list with `estimates`
# (theta hat), `beta`, `gamma` and `M`
sens_parts <- function(fit, data) {
  return(list(
    raw = reconstruct(fit, data, trace = FALSE),
    id = reconstruct(imr_identify(fit, data), data, trace = FALSE)
  ))
}

# Accuracy of a fit, as one row.
#   rrmse_all, rrmse_obs, rrmse_unobs: relative error of theta hat on all,
#     observed and unobserved entries (rmse_* likewise, not relative).
#   <C>_<estimate>_<target> for C in beta, gamma, M: relative error of the
#     raw or identified (id) estimate from the generating component (gen) or
#     from its part of T(theta) (T).
#   <C>_raw_id: how far the identification step moved the estimate.
#   rrmse_g01, ...: relative error of theta hat within each group of entries,
#     when `groups` (a matrix of group labels) is given.
sens_accuracy <- function(parts, d, groups = NULL) {
  err2 <- (parts$raw$estimates - d$theta)^2
  theta2 <- d$theta^2
  obs <- d$mask == 1
  row <- data.frame(
    rrmse_all = sqrt(sum(err2) / sum(theta2)),
    rrmse_obs = sqrt(sum(err2[obs]) / sum(theta2[obs])),
    rrmse_unobs = sqrt(sum(err2[!obs]) / sum(theta2[!obs])),
    rmse_all = sqrt(mean(err2)),
    rmse_obs = sqrt(mean(err2[obs])),
    rmse_unobs = sqrt(mean(err2[!obs]))
  )
  truth <- d$identified_truth
  for (name in c("beta", "gamma", "M")) {
    raw <- parts$raw[[name]]
    id <- parts$id[[name]]
    row[[paste0(name, "_raw_gen")]] <- rel_err(raw, d[[name]])
    row[[paste0(name, "_raw_T")]] <- rel_err(raw, truth[[name]])
    row[[paste0(name, "_id_gen")]] <- rel_err(id, d[[name]])
    row[[paste0(name, "_id_T")]] <- rel_err(id, truth[[name]])
    row[[paste0(name, "_raw_id")]] <- rel_err(raw, id)
  }
  if (!is.null(groups)) {
    by_group <- sqrt(tapply(err2, groups, sum) / tapply(theta2, groups, sum))
    row[sprintf("rrmse_g%02d", seq_along(by_group))] <- as.list(by_group)
  }
  return(row)
}

# Decile of every entry by its MAR score s_i + t_j: 1 = the entries least
# often observed under MAR, 10 = the most often. The observation probability
# of a MAR data set is increasing in the score, so its deciles are those of
# the score. The groups are the same for every mechanism of a replication.
sens_score_groups <- function(cfg, cell) {
  strengths <- cfg$missing$mar_strength
  cell$mechanism <- "mar"
  cell$mar_strength <- min(strengths[strengths > 0])
  prob <- sens_data(cfg, cell)$prob
  if (is.null(prob)) return(NULL)
  cuts <- quantile(prob, seq(0.1, 0.9, by = 0.1))
  return(matrix(findInterval(prob, cuts) + 1L, nrow(prob), ncol(prob)))
}

# The whole pipeline of one tuned cell: data, tuning, metrics.
# Returns list(row, history): one row of results, and the tuning path.
sens_cell <- function(cell, cfg, groups = NULL) {
  d <- sens_data(cfg, cell)
  if (!is.null(d$skipped)) {
    row <- cbind(cell, status = "skipped", message = d$skipped)
    return(list(row = row))
  }
  out <- sens_tune(d, cfg, cell$seed_fit)
  parts <- sens_parts(out$fit, out$data)
  row <- cbind(
    cell, status = "ok", sens_facts(d), sens_tuning(out),
    sens_accuracy(parts, d, groups)
  )
  history <- cbind(cell_id = cell$cell_id, out$history)
  return(list(row = row, history = history))
}
