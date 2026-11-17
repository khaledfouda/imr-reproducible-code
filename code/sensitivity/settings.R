# Settings of the sensitivity analysis. Every number of the study is here.
#
# size = "full"  : the design (about 22 to 35 core-hours).
# size = "draft" : fewer replications, enough to write the section. Its tasks
#                  are a subset of the full ones and share their folder, so a
#                  later "full" run only adds the missing tasks.
# size = "lite"  : a small problem that runs in minutes, to test the code.
sens_config <- function(size = c("draft", "full", "lite")) {
  size <- match.arg(size)
  cfg <- list(
    size = size,
    n_workers = 8,
    out_dir = "./output/sensitivity",
    # baseline data set: arguments of generate_simulated_data()
    base = list(n = 600, m = 600, r = 4, p = 3, q = 3, sparsity = 0.7, snr = 3),
    # tuned fits: arguments of imr(); the Lasso penalties are fixed at zero
    fit = list(
      center = FALSE, scale = TRUE, val_prop = 0.2, maxit = 2000, thresh = 1e-6
    ),
    # part 1: starting values and the monotone objective
    init = list(
      seeds = 1:10,                           # data sets
      n_random = 25,                          # random starts per data set
      thresh = c(tight = 1e-12, short = 1e-6),
      maxit = c(tight = 10000, short = 2000),
      rank_extra = 4,          # regime C: rank bound = tuned bound + 4
      lambda_factor = 1.1,     # regime C: lambda_m = 1.1^k * sigma_1(E)
      k_max = 40,
      cert_tol = 1e-6,         # certificate <= 1 + cert_tol
      increase_tol = 1e-12     # a relative increase of the objective counts
    ),
    # part 2: variability by source
    sources = list(
      seeds = 1:10, draws = 50, sources = c("split", "noise", "mask")
    ),
    # part 3: identification
    ident = list(reps = 500, n = c(200, 400, 600, 800, 1000)),
    # part 4: missingness mechanism (mar_strength = 0 is MCAR)
    missing = list(reps = 500, mar_strength = c(0, 0.25, 0.5, 1))
  )

  if (size == "draft") {
    cfg$init$seeds <- 1:5
    cfg$init$n_random <- 10
    cfg$sources$seeds <- 1:5
    cfg$ident$reps <- 9 * 3
    cfg$missing$reps <- 9 * 3
  }
  if (size == "lite") {
    cfg$n_workers <- 2
    cfg$out_dir <- "./output/sensitivity_lite"
    cfg$base$n <- 150
    cfg$base$m <- 150
    cfg$base$r <- 2
    cfg$init$seeds <- 1:2
    cfg$init$n_random <- 2
    cfg$init$maxit <- c(tight = 500, short = 200)
    cfg$sources$seeds <- 1
    cfg$sources$draws <- 3
    cfg$ident$reps <- 2
    cfg$ident$n <- c(100, 150)
    cfg$missing$reps <- 2
  }
  return(cfg)
}
