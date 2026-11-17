# Sensitivity analysis: the numbers behind each table and figure.
# Reads output/sensitivity/results_<part>.rds (written by 1_run.R) and saves
# output/sensitivity/summaries.rds. The tables and the figures themselves are
# made from these summaries.
library(dplyr)
library(tidyr)

out_dir <- "./output/sensitivity"
read_rows <- function(part) {
  readRDS(file.path(out_dir, paste0("results_", part, ".rds")))$rows
}

# Helpers ----
# Mean, standard deviation and Monte Carlo standard error (sd / sqrt(reps))
# of each variable, by group.
mc_summary <- function(rows, by, vars) {
  rows %>%
    pivot_longer(all_of(vars), names_to = "variable") %>%
    group_by(across(all_of(c(by, "variable")))) %>%
    summarise(
      mean = mean(value), sd = sd(value), mcse = sd / sqrt(n()), reps = n(),
      .groups = "drop"
    )
}

# Mean and standard error of the difference from the reference level `ref`
# of the column `by`, within each replication (`id`). Each replication needs
# exactly one row at the reference level.
paired_diff <- function(rows, id, by, ref, vars) {
  rows %>%
    pivot_longer(all_of(vars), names_to = "variable") %>%
    group_by(across(all_of(c(id, "variable")))) %>%
    mutate(diff = value - value[.data[[by]] == ref]) %>%
    group_by(across(all_of(c(by, "variable")))) %>%
    summarise(
      mean_diff = mean(diff), se_diff = sd(diff) / sqrt(n()), reps = n(),
      .groups = "drop"
    )
}

# Within-data-set standard deviation of each variable, pooled over the data
# sets (the square root of the mean within-data-set variance).
pooled_sd <- function(rows, vars) {
  rows %>%
    group_by(seed) %>%
    summarise(across(all_of(vars), var), .groups = "drop") %>%
    summarise(across(all_of(vars), ~ sqrt(mean(.x))))
}

# Part 4: missingness mechanism ----
summary_missing <- function(rows) {
  vars <- c("rrmse_all", "rrmse_unobs", "rrmse_obs", "beta_id_T",
            "gamma_id_T", "M_id_T", "lambda_m", "rank_fit")
  deciles <- grep("^rrmse_g", names(rows), value = TRUE)
  # the paired comparisons keep the replications with all mechanisms
  complete <- rows %>% group_by(task) %>% filter(all(status == "ok")) %>%
    ungroup()
  list(
    # how many draws were skipped (too few observed entries), by mechanism
    skipped = rows %>% group_by(mar_strength) %>%
      summarise(skipped = sum(status == "skipped"), reps = n()),
    # strength of the mechanism: pi_max / pi_min
    pi_ratio = complete %>% group_by(mar_strength) %>%
      summarise(pi_ratio = mean(pi_max / pi_min),
                min_row_obs = min(min_row_obs)),
    # table: mean, sd and Monte Carlo SE; figure, panel 1: rrmse_all
    table = mc_summary(complete, "mar_strength", vars),
    # table: paired differences from MCAR
    paired = paired_diff(complete, "task", "mar_strength", 0, vars),
    # figure, panel 2: error by decile of the MAR score
    deciles = mc_summary(complete, "mar_strength", deciles)
  )
}

# Part 3: identification ----
summary_ident <- function(rows) {
  rows <- rows %>% filter(status == "ok")
  errors <- function(name) {
    paste0(name, c("_raw_gen", "_raw_T", "_id_gen", "_id_T"))
  }
  list(
    # figure: error on beta against n, and the gap (dashed line)
    beta = mc_summary(rows, c("identified", "n"),
                      c(errors("beta"), "gap_beta")),
    # appendix table: gamma, M and the prediction error
    others = mc_summary(
      rows, c("identified", "n"),
      c(errors("gamma"), "gap_gamma", errors("M"), "gap_M", "rrmse_unobs")
    )
  )
}

# Part 2: variability by source ----
# One row per source: the pooled within-data-set standard deviation of each
# output. The "start" row comes from part 1 (regime T, short threshold) and
# the "everything" row from the MCAR replications of part 4.
summary_sources <- function(rows, init_rows, missing_rows) {
  rows <- rows %>% filter(status == "ok") %>%
    mutate(log_lambda_m = log(lambda_m))
  accuracy <- c("rrmse_unobs", "beta_id_T", "gamma_id_T", "M_id_T")
  tuning <- c("log_lambda_m", "rank_fit")

  by_source <- rows %>% group_by(source) %>%
    group_modify(~ pooled_sd(.x, c(accuracy, tuning))) %>% ungroup()
  start <- init_rows %>% filter(regime == "T", level == "short") %>%
    pooled_sd(accuracy) %>% mutate(source = "start")
  everything <- missing_rows %>%
    filter(status == "ok", mar_strength == 0) %>%
    mutate(log_lambda_m = log(lambda_m)) %>%
    summarise(across(all_of(c(accuracy, tuning)), sd)) %>%
    mutate(source = "everything")
  list(
    table = bind_rows(start, by_source, everything),
    # the grid is coarse: number of distinct values selected per data set
    distinct = rows %>% group_by(source, seed) %>%
      summarise(lambda_m = n_distinct(lambda_m),
                rank_fit = n_distinct(rank_fit), .groups = "drop")
  )
}

# Part 1: starting values and the monotone objective ----
summary_init <- function(rows) {
  # relative excess of the final objective over the best start of the same
  # data set, regime and threshold
  rows <- rows %>% group_by(seed, regime, level) %>%
    mutate(excess = (objective - min(objective)) / min(objective)) %>%
    ungroup()
  others <- rows %>% filter(start != "default")
  differences <- grep("^d_", names(rows), value = TRUE)
  list(
    # table: median and maximum over the data sets and the starts
    table = others %>%
      pivot_longer(all_of(c("excess", "certificate", differences)),
                   names_to = "variable") %>%
      group_by(regime, level, variable) %>%
      summarise(median = median(value), max = max(value), .groups = "drop"),
    # range of the test error across the starts of a data set
    test_range = rows %>% group_by(seed, regime, level) %>%
      summarise(range = max(rrmse_unobs) - min(rrmse_unobs),
                .groups = "drop") %>%
      group_by(regime, level) %>%
      summarise(median = median(range), max = max(range), .groups = "drop"),
    # iterations and convergence, by kind of start
    iterations = rows %>% group_by(regime, level, start) %>%
      summarise(median_iter = median(n_iter), max_iter = max(n_iter),
                share_at_cap = mean(!converged), .groups = "drop"),
    # monotone objective: the largest relative increase over all the fits,
    # and the number of iterations with an increase above the tolerance
    monotone = rows %>%
      summarise(max_rel_increase = max(max_rel_increase),
                n_increase = sum(n_increase), fits = n())
  )
}

# Objective traces of one data set, for the figure F_k - F_min against the
# iteration: one row per fit and iteration.
init_traces <- function(seed) {
  setup <- readRDS(file.path(out_dir, "init_setup",
                             sprintf("seed_%02d.rds", seed)))$extra
  # a data frame from the objective of one fit
  trace <- function(fit, regime, level, start) {
    data.frame(regime = regime, level = level, start = start,
               iteration = seq_along(fit$objective), objective = fit$objective)
  }
  default <- bind_rows(lapply(c("T", "C"), function(regime) {
    bind_rows(lapply(names(setup$fits[[regime]]), function(level) {
      trace(setup$fits[[regime]][[level]], regime, level, "default")
    }))
  }))
  files <- list.files(file.path(out_dir, "init"), full.names = TRUE,
                      pattern = sprintf("^seed_%02d_.*\\.rds$", seed))
  others <- bind_rows(lapply(files, function(file) {
    result <- readRDS(file)
    bind_rows(lapply(names(result$extra$fits), function(level) {
      row <- result$rows[result$rows$level == level, ]
      start <- paste(row$start, row$start_id)
      trace(result$extra$fits[[level]], row$regime, level, start)
    }))
  }))
  bind_rows(default, others)
}

# Run ----
missing_rows <- read_rows("missing")
ident_rows <- read_rows("ident")
sources_rows <- read_rows("sources")
init_rows <- bind_rows(read_rows("init_setup"), read_rows("init"))

summaries <- list(
  missing = summary_missing(missing_rows),
  ident = summary_ident(ident_rows),
  sources = summary_sources(sources_rows, init_rows, missing_rows),
  init = summary_init(init_rows),
  init_traces = init_traces(seed = 1)
)
saveRDS(summaries, file.path(out_dir, "summaries.rds"))
