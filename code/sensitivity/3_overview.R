# Sensitivity analysis: an overview of the results, to read them.
#
# Prints a few small tables per part and draws four figures. It reads
# output/sensitivity/results_<part>.rds (written by 1_run.R) and saves
# nothing. The figure of part 1 also needs summaries.rds (2_summaries.R).
# Run it from the root of the simulation repository, section by section in
# RStudio, or all at once with source("code/sensitivity/3_overview.R").
#
# All errors are relative: ||estimate - target||_F / ||target||_F.
library(dplyr)
library(tidyr)
library(ggplot2)

out_dir <- "./article/imr-reproducible-code/output/sensitivity"
read_results <- function(part) {
  readRDS(file.path(out_dir, paste0("results_", part, ".rds")))
}

# Helpers ----
# "mean (Monte Carlo standard error)" of a column, as text
mean_se <- function(x) {
  sprintf("%.4f (%.4f)", mean(x), sd(x) / sqrt(length(x)))
}

# print a table under its title
show <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(as.data.frame(table), digits = 3, row.names = FALSE)
}

# Part 4: missingness mechanism ----
# One replication = the same signal, noise and uniforms under four
# mechanisms: MCAR (mar_strength = 0) and MAR of increasing strength.
missing_rows <- read_results("missing")$rows %>% filter(status == "ok")

# pi_ratio = largest / smallest observation probability.
# all, unobserved, observed: error of the fitted matrix on those entries.
# beta, gamma, M: error of the identified estimates.
show(
  "Part 4. Errors by mechanism: mean (SE)",
  missing_rows %>% group_by(mar_strength) %>% summarise(
    reps = n(), pi_ratio = mean(pi_max / pi_min),
    fewest_in_a_row = min(min_row_obs),
    all = mean_se(rrmse_all), unobserved = mean_se(rrmse_unobs),
    observed = mean_se(rrmse_obs), beta = mean_se(beta_id_T),
    gamma = mean_se(gamma_id_T), M = mean_se(M_id_T)
  )
)

# The same errors minus their MCAR value in the same replication. The SE is
# small because the replications are paired.
paired <- missing_rows %>%
  group_by(task) %>% filter(any(mar_strength == 0)) %>%
  mutate(across(
    c(rrmse_all, rrmse_unobs, rrmse_obs, beta_id_T, gamma_id_T, M_id_T),
    ~ .x - .x[mar_strength == 0]
  )) %>%
  ungroup() %>% filter(mar_strength > 0)
show(
  "Part 4. Increase over MCAR, within a replication: mean (SE)",
  paired %>% group_by(mar_strength) %>% summarise(
    all = mean_se(rrmse_all), unobserved = mean_se(rrmse_unobs),
    observed = mean_se(rrmse_obs), beta = mean_se(beta_id_T),
    gamma = mean_se(gamma_id_T), M = mean_se(M_id_T)
  )
)

# The entries are cut in ten groups by their MAR score: group 1 is the
# least often observed under MAR, group 10 the most often.
deciles <- missing_rows %>%
  pivot_longer(starts_with("rrmse_g"), names_to = "decile",
    names_prefix = "rrmse_g", values_to = "error") %>%
  group_by(mar_strength, decile = as.integer(decile)) %>%
  summarise(error = mean(error), .groups = "drop")
show(
  "Part 4. Error by decile of the MAR score (1 = least often observed)",
  deciles %>% pivot_wider(names_from = decile, values_from = error)
)
print(
  ggplot(deciles, aes(decile, error, colour = factor(mar_strength))) +
    geom_line() + geom_point() +
    scale_x_continuous(breaks = 1:10) +
    labs(title = "Part 4: error by decile of the MAR score",
      x = "decile (1 = least often observed)",
      y = "relative error of the fitted matrix",
      colour = "MAR strength")
)

# Part 3: identification ----
# One replication = five values of n, each with the truth as generated (it
# violates the conditions) and with the truth projected to satisfy them.
ident_rows <- read_results("ident")$rows %>% filter(status == "ok") %>%
  mutate(truth = ifelse(identified, "identified", "as generated"))

# gap                  : distance from the generating beta to the beta of
#                        the identified truth (0 when the truth is
#                        identified). No estimator of the identified beta
#                        can be closer than this to the generating beta.
# raw_vs_generating    : the raw estimate against the generating beta.
# id_vs_generating     : the identified estimate against the generating beta.
# id_vs_identified     : the identified estimate against its own target.
show(
  "Part 3. Error on beta, by truth and n: mean",
  ident_rows %>% group_by(truth, n) %>% summarise(
    reps = n(), gap = mean(gap_beta),
    raw_vs_generating = mean(beta_raw_gen),
    id_vs_generating = mean(beta_id_gen),
    id_vs_identified = mean(beta_id_T),
    se = sd(beta_id_T) / sqrt(n()), .groups = "drop"
  )
)
show(
  "Part 3. Other identified components and the test error: mean",
  ident_rows %>% group_by(truth, n) %>% summarise(
    gamma = mean(gamma_id_T), M = mean(M_id_T),
    test_error = mean(rrmse_unobs), .groups = "drop"
  )
)

# The rank bound selected by the tuning and the rank of the fit. When the
# two are always equal, the fitted rank is the bound and not an estimate.
show(
  "Part 3. Selected rank bound and fitted rank: number of fits",
  ident_rows %>% count(n, rank_bound = rank_in, rank_fit, name = "fits")
)

beta_lines <- ident_rows %>%
  pivot_longer(c(beta_raw_gen, beta_id_gen, beta_id_T, gap_beta),
    names_to = "line", values_to = "error") %>%
  group_by(truth, n, line) %>%
  summarise(error = mean(error), .groups = "drop") %>%
  mutate(line = recode(
    line,
    beta_raw_gen = "raw estimate, generating beta",
    beta_id_gen = "identified estimate, generating beta",
    beta_id_T = "identified estimate, identified beta",
    gap_beta = "gap between the two targets"
  ))
print(
  ggplot(beta_lines, aes(n, error, colour = line)) +
    geom_line() + geom_point() +
    facet_wrap(~truth, labeller = label_both) +
    labs(title = "Part 3: error on beta against n",
      y = "relative error", colour = "estimate, target") +
    theme(legend.position = "bottom", legend.direction = "vertical")
)

# Part 1: starting values and the monotone objective ----
# Read before part 2, which uses its rows.
# Regime T: the tuned lambda_m and rank bound. Regime C: a larger lambda_m
# and bound, where the fit minimizes the convex problem (certificate = 1).
# Starts: default, zero, and random. Thresholds: tight and short.
setup <- read_results("init_setup")
init_rows <- bind_rows(setup$rows, read_results("init")$rows) %>%
  group_by(seed, regime, level) %>%
  mutate(
    # final objective, relative to the best start of the same data set
    excess = (objective - min(objective)) / min(objective),
    # range of the test error across the starts of the same data set
    test_range = max(rrmse_unobs) - min(rrmse_unobs)
  ) %>%
  ungroup()

# found = FALSE would mean that the search of regime C failed
show(
  "Part 1. The two regimes of each data set",
  bind_rows(lapply(setup$extra, `[[`, "regimes"), .id = "data_set")
)

# certificate : sigma_1(residual) / lambda_m; 1 = minimizer of the convex
#               problem, above 1 = the rank bound is binding.
# at_cap      : fits stopped by the iteration cap.
# excess      : largest relative excess of the final objective.
show(
  "Part 1. Fits by regime, threshold and start",
  init_rows %>% group_by(regime, level, start) %>% summarise(
    fits = n(), iterations = median(n_iter), at_cap = sum(!converged),
    certificate = median(certificate),
    rank_fit = paste(unique(range(rank_fit)), collapse = " to "),
    rank_bound = max(rank), excess = max(excess),
    test_error = mean(rrmse_unobs), .groups = "drop"
  )
)

# Relative difference between the fit of a zero or random start and the
# default-start fit of the same data set, regime and threshold.
# far : starts whose fitted matrix differs by more than 0.1%.
show(
  "Part 1. Difference from the default-start fit (zero and random starts)",
  init_rows %>% filter(start != "default") %>%
    group_by(regime, level) %>% summarise(
      fits = n(), theta_median = median(d_theta_all),
      theta_max = max(d_theta_all), far = sum(d_theta_all > 1e-3),
      beta_id_max = max(d_beta_id), beta_raw_max = max(d_beta_raw),
      M_id_max = max(d_M_id), test_range_max = max(test_range),
      .groups = "drop"
    )
)
show(
  "Part 1. Tight threshold: starts far from the default-start fit, by data set",
  init_rows %>% filter(start != "default", level == "tight") %>%
    group_by(regime, seed) %>%
    summarise(far = sum(d_theta_all > 1e-3), starts = n(),
      theta_max = max(d_theta_all), excess_max = max(excess),
      .groups = "drop")
)
cat(sprintf(
  "\nMonotone objective: %d fits, %d iterations with an increase, largest relative change %.1e (negative = decrease).\n",
  nrow(init_rows), sum(init_rows$n_increase), max(init_rows$max_rel_increase)
))

# Objective minus the smallest value reached in the same regime and
# threshold, against the iteration, for the first data set (log scales).
traces <- readRDS(file.path(out_dir, "summaries.rds"))$init_traces %>%
  group_by(regime, level) %>%
  mutate(gap = objective - min(objective),
    kind = sub(" .*", "", start)) %>%
  ungroup() %>% filter(gap > 0)
print(
  ggplot(traces, aes(iteration, gap, group = start, colour = kind)) +
    geom_line(alpha = 0.6) +
    scale_x_log10() + scale_y_log10() +
    facet_grid(regime ~ level, labeller = label_both) +
    labs(title = "Part 1: objective minus its smallest value, data set 1",
      y = "objective - minimum", colour = "start")
)

# Part 2: variability by source ----
# From one data set, one thing is redrawn 50 times: the train/validation
# split, the noise, or the mask. The table gives the standard deviation of
# each output within a data set, pooled over the data sets.
sources_rows <- read_results("sources")$rows %>% filter(status == "ok")
outputs <- c("rrmse_unobs", "beta_id_T", "gamma_id_T", "M_id_T")

pooled_sd <- function(rows) {
  rows %>% group_by(seed) %>%
    summarise(across(all_of(outputs), var), .groups = "drop") %>%
    summarise(across(all_of(outputs), ~ sqrt(mean(.x))))
}
mcar_rows <- missing_rows %>% filter(mar_strength == 0)
# rows: the three sources; the start (part 1, regime T, short threshold);
# everything redrawn (the MCAR replications of part 4); the mean, for scale
show(
  "Part 2. Standard deviation of each output, by what is redrawn",
  bind_rows(
    sources_rows %>% group_by(source) %>% group_modify(~ pooled_sd(.x)) %>%
      ungroup(),
    init_rows %>% filter(regime == "T", level == "short") %>% pooled_sd() %>%
      mutate(source = "start"),
    mcar_rows %>% summarise(across(all_of(outputs), sd)) %>%
      mutate(source = "everything"),
    mcar_rows %>% summarise(across(all_of(outputs), mean)) %>%
      mutate(source = "(mean, for scale)")
  )
)

# what the tuning selected across the draws: lambda_sd is the standard
# deviation of lambda_m within a data set, pooled over the data sets
show(
  "Part 2. Tuning across the draws: lambda_m and the fitted ranks",
  sources_rows %>% group_by(source) %>% summarise(
    lambda_mean = mean(lambda_m),
    lambda_sd = sqrt(mean(tapply(lambda_m, seed, var))),
    ranks_fitted = paste(sort(unique(rank_fit)), collapse = ", ")
  )
)
print(
  sources_rows %>% group_by(seed, source) %>%
    mutate(deviation = rrmse_unobs - mean(rrmse_unobs)) %>%
    ggplot(aes(source, deviation)) +
    geom_boxplot() +
    labs(title = "Part 2: test error minus its mean in the data set",
      x = "what is redrawn", y = "deviation of the test error")
)
