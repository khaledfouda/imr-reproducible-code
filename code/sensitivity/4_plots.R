# Sensitivity analysis: the four figures of the paper, one per case.
#
# Reads output/sensitivity/results_<part>.rds (written by 1_run.R) and writes
# the figures (.pdf) to output/sensitivity/paper/. It ends by printing the
# numbers quoted in the text of the section.
# Run from the root of the simulation repository:
#   Rscript code/sensitivity/4_plots.R
# The theme and the colours are those of code/simulation/2_tables_plots.R.
library(dplyr)
library(tidyr)
library(ggplot2)
library(RSSthemes)

out_dir <- "./output/sensitivity"
paper_dir <- file.path(out_dir, "paper")
dir.create(paper_dir, showWarnings = FALSE)
read_results <- function(part) {
  readRDS(file.path(out_dir, paste0("results_", part, ".rds")))
}

# Helpers ----
# "mean (standard deviation)" of a column, as text
mean_sd <- function(x, digits = 3) {
  sprintf("%.*f (%.*f)", digits, mean(x), digits, sd(x))
}

# "mean (Monte Carlo standard error)" of a column, as text
mean_se <- function(x, digits = 4) {
  sprintf("%.*f (%.*f)", digits, mean(x), digits, sd(x) / sqrt(length(x)))
}

save_figure <- function(plot, name, width = 7, height = 4.5) {
  ggsave(file.path(paper_dir, paste0(name, ".pdf")), plot,
         device = cairo_pdf, width = width, height = height, units = "in")
}

theme_paper <- theme_significance() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold", size = 11))

# Three colours that colour-blind readers can tell apart, and four steps of
# one hue for the ordered strengths of the missingness mechanism.
series_colors <- c(signif_blue, signif_orange, signif_red)
strength_colors <- c("#4fbcd8", signif_blue, "#00718e", "#00485b")

# Missingness mechanism ----
# The comparisons are paired: only the replicates with all four mechanisms.
missing_rows <- read_results("missing")$rows %>%
  group_by(task) %>% filter(all(status == "ok")) %>% ungroup() %>%
  mutate(ratio = pi_max / pi_min)

mechanisms <- missing_rows %>% group_by(mar_strength) %>%
  summarise(ratio = mean(ratio), pi_min = mean(pi_min), pi_max = mean(pi_max),
            lambda_m = mean(lambda_m), fewest_in_a_row = min(min_row_obs),
            replicates = n(), .groups = "drop") %>%
  mutate(label = ifelse(mar_strength == 0, "MCAR (1)",
                        sprintf("MAR (%.0f)", ratio)))
missing_rows <- missing_rows %>%
  mutate(mechanism = factor(mar_strength, levels = mechanisms$mar_strength,
                            labels = mechanisms$label))

# Figure: error of the fitted matrix within each decile of the score that
# drives the MAR mechanisms, one line per mechanism
by_decile <- missing_rows %>%
  pivot_longer(starts_with("rrmse_g"), names_to = "decile",
               names_prefix = "rrmse_g", values_to = "rrmse") %>%
  group_by(mechanism, decile = as.integer(decile)) %>%
  summarise(mean = mean(rrmse), sd = sd(rrmse), .groups = "drop")
missing_figure <- ggplot(by_decile, aes(decile, mean, color = mechanism,
                                        shape = mechanism)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = strength_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Decile of the observation score (1 = least often observed)",
       y = expression(RRMSE(Theta)),
       color = expression("Mechanism" ~ (pi[max] / pi[min]) * ":"),
       shape = expression("Mechanism" ~ (pi[max] / pi[min]) * ":")) +
  theme_paper
save_figure(missing_figure, "figure_sens_missing")

# Identified components ----
ident_rows <- read_results("ident")$rows %>% filter(status == "ok") %>%
  mutate(truth = factor(
    identified, levels = c(TRUE, FALSE),
    labels = c("Truth satisfies (C1)–(C3)",
               "Truth violates (C1)–(C3)")
  ))

# Figure: error on beta against n. With a truth that satisfies the
# conditions the two targets are the same, so one line is enough.
line_labels <- c(
  beta_id_T = "Identified estimate vs. identified truth",
  beta_id_gen = "Identified estimate vs. generating coefficients",
  beta_raw_gen = "Raw estimate vs. generating coefficients"
)
beta_lines <- ident_rows %>%
  pivot_longer(all_of(names(line_labels)), names_to = "line",
               values_to = "error") %>%
  filter(!(identified & line == "beta_id_gen")) %>%
  group_by(truth, n, line) %>%
  summarise(mean = mean(error), sd = sd(error), .groups = "drop") %>%
  mutate(line = factor(line, levels = names(line_labels),
                       labels = line_labels))
ident_figure <- ggplot(beta_lines, aes(n, mean, color = line, fill = line,
                                       shape = line)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15,
              color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = series_colors) +
  scale_fill_manual(values = series_colors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_x_continuous(breaks = sort(unique(beta_lines$n))) +
  facet_wrap(~truth) +
  labs(x = "Number of rows (n)", y = expression(RRMSE(beta)), color = NULL,
       fill = NULL, shape = NULL) +
  theme_paper
save_figure(ident_figure, "figure_sens_identification", width = 12,
            height = 5)

# Initialization ----
# Only the fits at the tuned parameters (regime "T" of init.R) are used.
# Starts: default, zero and random. Tolerances: tight (1e-12) and short
# (1e-6, the tolerance of the other cases).
init_rows <- bind_rows(read_results("init_setup")$rows,
                       read_results("init")$rows) %>%
  filter(regime == "T") %>%
  group_by(seed, level) %>%
  mutate(
    # final objective, relative to the best start of the same data set
    excess = (objective - min(objective)) / min(objective),
    # range of the test error across the starts of the same data set
    test_range = max(rrmse_unobs) - min(rrmse_unobs)
  ) %>%
  ungroup() %>%
  mutate(kind = factor(start, levels = c("default", "zero", "random"),
                       labels = c("Default", "Zero", "Random")))

# Figure: the objective after each iteration for one data set (tight
# tolerance), minus the smallest value reached from any start
trace_seed <- 1
traces_of <- function(seed) {
  name <- sprintf("seed_%02d", seed)
  default <- readRDS(
    file.path(out_dir, "init_setup", paste0(name, ".rds"))
  )$extra$fits$T$tight
  trace <- function(fit, start) {
    data.frame(start = start, iteration = seq_along(fit$objective),
               objective = fit$objective)
  }
  files <- list.files(file.path(out_dir, "init"), full.names = TRUE,
                      pattern = paste0("^", name, "_T_"))
  others <- lapply(files, function(file) {
    result <- readRDS(file)
    row <- result$rows[1, ]
    trace(result$extra$fits$tight, paste(row$start, row$start_id))
  })
  bind_rows(trace(default, "default"), others)
}
traces <- traces_of(trace_seed) %>%
  mutate(gap = objective - min(objective),
         kind = factor(sub(" .*", "", start),
                       levels = c("default", "zero", "random"),
                       labels = c("Default", "Zero", "Random"))) %>%
  filter(gap > 0)
init_figure <- ggplot(traces, aes(iteration, gap, group = start,
                                  color = kind)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  scale_color_manual(values = series_colors) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_log()) +
  labs(x = "Iteration (log scale)",
       y = "Objective minus its smallest value (log scale)",
       color = "Initialization:") +
  theme_paper
save_figure(init_figure, "figure_sens_init")

# Sources of variability ----
# From one data set, one thing is redrawn: the starting values (the fits
# above, at the short tolerance), the validation split, the noise or the
# set of observed entries. The figure shows the test error minus its mean
# in the data set.
source_labels <- c(start = "Initial values", split = "Validation split",
                   noise = "Noise", mask = "Observed entries")
sources_rows <- bind_rows(
  read_results("sources")$rows %>% filter(status == "ok"),
  init_rows %>% filter(level == "short") %>% mutate(source = "start")
) %>%
  mutate(source = factor(source, levels = names(source_labels),
                         labels = source_labels)) %>%
  group_by(seed, source) %>%
  mutate(deviation = rrmse_unobs - mean(rrmse_unobs)) %>%
  ungroup()
sources_figure <- ggplot(sources_rows, aes(source, deviation)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_boxplot(color = signif_blue, fill = signif_blue, alpha = 0.15,
               outlier.size = 1) +
  labs(x = "Component redrawn",
       y = "Test RRMSE minus its mean in the data set") +
  theme_paper
save_figure(sources_figure, "figure_sens_sources")

# Numbers quoted in the text ----
errors <- c("rrmse_all", "rrmse_unobs", "rrmse_obs", "beta_id_T",
            "gamma_id_T", "M_id_T")
cat("\n== Missingness: the mechanisms ==\n")
print(as.data.frame(mechanisms), digits = 3)
cat("errors: mean (sd)\n")
print(as.data.frame(missing_rows %>% group_by(mar_strength) %>%
  summarise(across(all_of(errors), mean_sd))), row.names = FALSE)
cat("difference from MCAR within a replicate: mean (standard error)\n")
print(as.data.frame(missing_rows %>% group_by(task) %>%
  mutate(across(all_of(errors), ~ .x - .x[mar_strength == 0])) %>%
  ungroup() %>% filter(mar_strength > 0) %>% group_by(mar_strength) %>%
  summarise(across(all_of(errors), mean_se))), row.names = FALSE)
cat("first and last deciles: mean\n")
print(as.data.frame(by_decile %>% filter(decile %in% c(1, 10)) %>%
  select(-sd) %>% pivot_wider(names_from = decile, values_from = mean,
                              names_prefix = "decile_")),
  digits = 3, row.names = FALSE)

cat("\n== Identified components: mean errors ==\n")
print(as.data.frame(ident_rows %>% group_by(truth, n) %>% summarise(
  replicates = n(), beta_raw_gen = mean(beta_raw_gen),
  beta_id_gen = mean(beta_id_gen), beta_id_T = mean(beta_id_T),
  targets_apart = mean(gap_beta), gamma_id_T = mean(gamma_id_T),
  M_id_T = mean(M_id_T), test = mean(rrmse_unobs), .groups = "drop"
)), digits = 3, row.names = FALSE)

cat("\n== Initialization (tuned parameters) ==\n")
cat("data sets:", n_distinct(init_rows$seed), " random starts per data set:",
    n_distinct(init_rows$start_id[init_rows$start == "random"]), "\n")
cat("iterations by tolerance and start (median), fits at the cap:\n")
print(as.data.frame(init_rows %>% group_by(level, kind) %>% summarise(
  fits = n(), iterations = median(n_iter), at_cap = sum(!converged),
  .groups = "drop"
)), row.names = FALSE)
cat("tight tolerance, by data set: differences from the default start\n")
print(as.data.frame(init_rows %>% filter(level == "tight") %>%
  group_by(seed) %>% summarise(
    far = sum(d_theta_all > 1e-3, na.rm = TRUE),
    theta_max = max(d_theta_all, na.rm = TRUE),
    beta_raw_max = max(d_beta_raw, na.rm = TRUE),
    beta_id_max = max(d_beta_id, na.rm = TRUE),
    objective_max = max(excess), test_min = min(rrmse_unobs),
    test_max = max(rrmse_unobs), .groups = "drop"
  )), digits = 3, row.names = FALSE)
cat("short tolerance: differences from the default start\n")
print(as.data.frame(init_rows %>% filter(level == "short") %>% summarise(
  theta_median = median(d_theta_all, na.rm = TRUE),
  theta_max = max(d_theta_all, na.rm = TRUE), test_range = max(test_range)
)), digits = 3, row.names = FALSE)
cat(sprintf(
  "monotone objective: %d fits, %d iterations, %d with an increase, largest relative change %.1e\n",
  nrow(init_rows), sum(init_rows$n_iter), sum(init_rows$n_increase),
  max(init_rows$max_rel_increase)
))

cat("\n== Sources: standard deviation within a data set, pooled ==\n")
print(as.data.frame(sources_rows %>% group_by(source, seed) %>%
  summarise(variance = var(rrmse_unobs), draws = n(), .groups = "drop") %>%
  group_by(source) %>%
  summarise(sd = sqrt(mean(variance)), data_sets = n(),
            draws = first(draws))), digits = 3, row.names = FALSE)
mcar_rows <- missing_rows %>% filter(mar_strength == 0)
cat(sprintf(
  "all components redrawn (MCAR replicates): mean %.4f, sd %.4f; lambda_m sd / mean within a data set: %.4f\n",
  mean(mcar_rows$rrmse_unobs), sd(mcar_rows$rrmse_unobs),
  max(tapply(read_results("sources")$rows$lambda_m,
             read_results("sources")$rows$task, sd) /
        tapply(read_results("sources")$rows$lambda_m,
               read_results("sources")$rows$task, mean))
))
