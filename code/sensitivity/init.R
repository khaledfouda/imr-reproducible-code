# Part 1: starting values, and the monotone objective.
#
# For each data set, the tuning parameters are fixed in two regimes:
#   T (tuned) : lambda_m and the rank bound selected by the tuning;
#   C (convex): a larger rank bound, and the smallest lambda_m on a grid for
#               which the certificate sigma_1(E) / lambda_m is 1 and the
#               fitted rank is below the bound. The fit is then a minimizer
#               of the convex problem.
# Every start is run at two thresholds ("tight" and "short").
#
# Two stages:
#   init_setup_block(): one task per data set. Tuning, the two regimes, and
#     the fits from the default start (the reference of all comparisons).
#   init_block(): one task per data set, regime and start. The fits from the
#     zero start and from the random starts, compared with the reference.

# Cells of the second stage. start_id = 0 is the zero start (zero
# coefficients, truncated SVD of the zero-filled Y); 1, 2, ... are random.
design_init <- function(cfg) {
  cells <- expand.grid(
    start_id = 0:cfg$init$n_random, regime = c("T", "C"),
    seed = cfg$init$seeds, stringsAsFactors = FALSE
  )
  cells$start <- ifelse(cells$start_id == 0, "zero", "random")
  cells$seed_start <- 10000 * cells$seed + cells$start_id
  cells$task <- sprintf(
    "seed_%02d_%s_%s_%02d", cells$seed, cells$regime, cells$start,
    cells$start_id
  )
  return(cells)
}

# Certificate sigma_1(E) / lambda_m, with E the residual at the observed
# entries. At most 1: the fit minimizes the convex problem.
init_certificate <- function(fit, lambda_m) {
  return(svd_opt(fit$residuals, 2)$d[1] / lambda_m)
}

# Random start with exactly `rank` columns: U and V are orthonormal bases of
# random subspaces, d holds the leading singular values of the zero-filled Y,
# and the coefficients are N(0, 1) on the scale of the stored designs.
# The two deterministic starts need no warm start (NULL).
init_start <- function(cell, data, d_scale) {
  if (cell$start != "random") return(NULL)
  n <- nrow(data$Y)
  m <- ncol(data$Y)
  p <- ncol(data$X)
  q <- ncol(data$Z)
  r <- cell$rank
  set.seed(cell$seed_start)
  return(list(
    u = qr.Q(qr(matrix(rnorm(n * r), n, r))),
    d = d_scale[seq_len(r)],
    v = qr.Q(qr(matrix(rnorm(m * r), m, r))),
    beta = matrix(rnorm(p * m), p, m),
    gamma = matrix(rnorm(n * q), n, q),
    beta0 = NULL, gamma0 = NULL
  ))
}

# Differences between a fit and the reference fit, as one row: theta hat on
# all, observed and unobserved entries, then the raw and the identified
# components.
init_differences <- function(parts, ref, obs) {
  est <- parts$raw$estimates
  ref_est <- ref$raw$estimates
  row <- data.frame(
    d_theta_all = rel_err(est, ref_est),
    d_theta_obs = rel_err(est[obs], ref_est[obs]),
    d_theta_unobs = rel_err(est[!obs], ref_est[!obs])
  )
  for (name in c("beta", "gamma", "M")) {
    row[[paste0("d_", name, "_raw")]] <-
      rel_err(parts$raw[[name]], ref$raw[[name]])
    row[[paste0("d_", name, "_id")]] <-
      rel_err(parts$id[[name]], ref$id[[name]])
  }
  return(row)
}

# One fit at fixed tuning parameters, and its row of results.
#   cell      : seed, regime, rank, lambda_m, start, start_id, seed_start.
#   level     : "tight" or "short" (threshold and iteration cap).
#   ref_parts : sens_parts() of the reference fit, or NULL for the reference.
# Returns list(row, fit). The fit holds the objective after each iteration
# (fit$objective) and the coefficients.
init_fit_row <- function(data, d, cell, warm_start, level, ref_parts, cfg) {
  convergence <- imr_convergence(
    maxit = cfg$init$maxit[[level]], thresh = cfg$init$thresh[[level]],
    ls_initial = cell$start != "zero"
  )
  # the default starts draw random numbers (irlba), so the fit is seeded
  set.seed(cell$seed_start)
  start_time <- Sys.time()
  fit <- suppressMessages(imr_fit(
    data, rank = cell$rank, lambda_m = cell$lambda_m,
    convergence = convergence, warm_start = warm_start,
    return_objective = TRUE
  ))
  secs <- as.numeric(Sys.time() - start_time, units = "secs")

  objective <- fit$objective
  # relative change of the objective at each iteration: positive = increase
  steps <- diff(objective) / objective[-length(objective)]
  parts <- sens_parts(fit, data)
  row <- cbind(
    cell, level = level,
    thresh = cfg$init$thresh[[level]], maxit = cfg$init$maxit[[level]],
    objective = objective[length(objective)],
    max_rel_increase = max(steps),
    n_increase = sum(steps > cfg$init$increase_tol),
    certificate = init_certificate(fit, cell$lambda_m),
    rank_fit = sum(fit$coefficients$d > 0),
    n_iter = fit$meta$n_iter, converged = fit$meta$converged, time = secs,
    sens_accuracy(parts, d)
  )
  if (!is.null(ref_parts)) {
    row <- cbind(row, init_differences(parts, ref_parts, d$mask == 1))
  }
  return(list(row = row, fit = fit))
}

# Regime C: the smallest lambda_m = factor^k * sigma_1(E_T), k = 1, 2, ...,
# whose default-start fit has certificate 1 and a rank below the bound. E_T
# is the residual of the tuned fit. Returns lambda_m and the search path.
init_regime_c <- function(data, fit_t, rank_c, seed, cfg) {
  sigma_1 <- init_certificate(fit_t, 1)
  convergence <- imr_convergence(
    maxit = cfg$init$maxit[["tight"]], thresh = cfg$init$thresh[["tight"]]
  )
  path <- NULL
  for (k in seq_len(cfg$init$k_max)) {
    lambda_m <- cfg$init$lambda_factor^k * sigma_1
    set.seed(seed)
    fit <- suppressMessages(imr_fit(
      data, rank = rank_c, lambda_m = lambda_m, convergence = convergence
    ))
    step <- data.frame(
      k = k, lambda_m = lambda_m,
      certificate = init_certificate(fit, lambda_m),
      rank_fit = sum(fit$coefficients$d > 0), n_iter = fit$meta$n_iter
    )
    path <- rbind(path, step)
    found <- step$certificate <= 1 + cfg$init$cert_tol &&
      step$rank_fit < rank_c
    if (found) break
  }
  return(list(lambda_m = lambda_m, found = found, path = path))
}

# Stage 1 worker: one data set (task = its seed).
init_setup_block <- function(seed, cfg) {
  d <- sens_data(cfg, data.frame(seed = seed))
  out <- sens_tune(d, cfg, seed)
  data <- out$data
  levels <- names(cfg$init$thresh)

  # default-start fits of one regime, at both thresholds
  reference <- function(regime, rank, lambda_m) {
    cell <- data.frame(
      seed = seed, regime = regime, rank = rank, lambda_m = lambda_m,
      start = "default", start_id = 0, seed_start = 10000 * seed
    )
    fits <- lapply(levels, function(level) {
      init_fit_row(data, d, cell, NULL, level, NULL, cfg)
    })
    names(fits) <- levels
    return(fits)
  }

  # regime T: the tuned values
  rank_t <- out$params$rank_in
  ref_t <- reference("T", rank_t, out$params$lambda_m)
  # regime C: found from the residual of the tight fit of regime T
  rank_c <- rank_t + cfg$init$rank_extra
  search <- init_regime_c(data, ref_t$tight$fit, rank_c, seed, cfg)
  ref_c <- reference("C", rank_c, search$lambda_m)

  regimes <- data.frame(
    regime = c("T", "C"), rank = c(rank_t, rank_c),
    lambda_m = c(out$params$lambda_m, search$lambda_m),
    found = c(TRUE, search$found)
  )
  # scale of the random starts: leading singular values of the zero-filled Y
  set.seed(seed)
  d_scale <- svd_opt(data$Y, rank_c)$d

  return(list(
    rows = dplyr::bind_rows(
      lapply(c(ref_t, ref_c), `[[`, "row")
    ),
    history = cbind(seed = seed, out$history),
    extra = list(
      data = data, regimes = regimes, d_scale = d_scale, search = search$path,
      # the reference fits, complete: the second stage compares with them
      fits = list(
        T = lapply(ref_t, `[[`, "fit"), C = lapply(ref_c, `[[`, "fit")
      )
    )
  ))
}

# Stage 2 worker: one start of one data set in one regime (task = its cell),
# at both thresholds.
init_block <- function(cell, cfg, setup_dir) {
  setup_file <- file.path(setup_dir, sprintf("seed_%02d.rds", cell$seed))
  setup <- readRDS(setup_file)$extra
  data <- setup$data
  d <- sens_data(cfg, cell)
  regime <- setup$regimes[setup$regimes$regime == cell$regime, ]
  cell$rank <- regime$rank
  cell$lambda_m <- regime$lambda_m
  warm_start <- init_start(cell, data, setup$d_scale)

  levels <- names(cfg$init$thresh)
  results <- lapply(levels, function(level) {
    ref_parts <- sens_parts(setup$fits[[cell$regime]][[level]], data)
    result <- init_fit_row(data, d, cell, warm_start, level, ref_parts, cfg)
    # the residuals are dropped to keep the file small: reconstruct() still
    # gives theta hat and the components from the saved fit
    result$fit$residuals <- NULL
    return(result)
  })
  names(results) <- levels
  return(list(
    rows = dplyr::bind_rows(lapply(results, `[[`, "row")),
    extra = list(fits = lapply(results, `[[`, "fit"))
  ))
}
