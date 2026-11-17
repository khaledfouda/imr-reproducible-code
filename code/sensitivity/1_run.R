# Sensitivity analysis: run the simulations and save the results.
#
# Run from the root of the simulation repository, in a terminal, with one
# BLAS thread per worker:
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/1_run.R lite
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/1_run.R draft
#   OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/1_run.R full missing ident
# First argument: the size (lite, draft or full; see settings.R).
# Other arguments: the parts to run (default: all four).
# A run that is stopped and started again continues where it stopped.

args <- commandArgs(trailingOnly = TRUE)
size <- if (length(args) > 0) args[1] else "draft"
parts <- if (length(args) > 1) args[-1] else
  c("missing", "ident", "sources", "init")

# Setup ----
# The development version of IMR (the installed one is 1.0.0). Replace this
# line by library(IMR) once the new version is installed.
# devtools::load_all("../..", quiet = TRUE)
library(IMR)

# generate_simulated_data() is defined in code/simulation/helper.R, which
# also attaches the plotting packages. Only the function is evaluated here.
for (e in parse("./code/simulation/helper.R")) {
  is_generator <- is.call(e) && identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name("generate_simulated_data"))
  if (is_generator) eval(e)
}

source("./code/sensitivity/settings.R")
source("./code/sensitivity/runner.R")
source("./code/sensitivity/helpers.R")
source("./code/sensitivity/tuned.R")
source("./code/sensitivity/init.R")

cfg <- sens_config(size)
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  warning("OPENBLAS_NUM_THREADS is not 1: each worker may use every core.")
}

# Run the tasks of one part in parallel, then save everything in one file:
# rows (one per fit), history (the tuning paths), extra (per task), and the
# settings of the run.
run_part <- function(part, tasks, worker, ...) {
  dir <- file.path(cfg$out_dir, part)
  run_tasks(tasks, worker, dir, cfg$n_workers, cfg = cfg, ...)
  results <- collect_tasks(dir)
  results$config <- cfg
  results$info <- list(
    date = Sys.time(), imr_version = as.character(packageVersion("IMR")),
    r_version = R.version.string,
    n_tasks = length(tasks), n_done = length(results$extra)
  )
  saveRDS(results, file.path(cfg$out_dir, paste0("results_", part, ".rds")))
}

# Part 4: missingness mechanism ----
if ("missing" %in% parts) {
  cells <- design_missing(cfg)
  run_part("missing", split(cells, cells$task), tuned_block, deciles = TRUE)
}

# Part 3: identification ----
if ("ident" %in% parts) {
  cells <- design_ident(cfg)
  run_part("ident", split(cells, cells$task), tuned_block)
}

# Part 2: variability by source ----
if ("sources" %in% parts) {
  cells <- design_sources(cfg)
  run_part("sources", split(cells, cells$task), tuned_block)
}

# Part 1: starting values and the monotone objective ----
if ("init" %in% parts) {
  seeds <- as.list(cfg$init$seeds)
  names(seeds) <- sprintf("seed_%02d", cfg$init$seeds)
  run_part("init_setup", seeds, init_setup_block)

  cells <- design_init(cfg)
  run_part("init", split(cells, cells$task), init_block,
    setup_dir = file.path(cfg$out_dir, "init_setup"))
}
