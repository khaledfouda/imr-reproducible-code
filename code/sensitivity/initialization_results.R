# Sensitivity analysis, initialization: tables and figures.
#
# Reads results_initialization.rds (written by initialization_run.R), prints
# the tables and saves the figures in the folder paper/ of the run.
# Run from the root of the simulation repository:
#   Rscript code/sensitivity/initialization_results.R         (draft or full)
#   Rscript code/sensitivity/initialization_results.R lite    (the test run)
#
# One row of the results is one fit: one start of one version of a signal at
# one threshold. Nothing is tuned, so the fits of a version of a signal
# differ in their start only. An error is the RRMSE of the paper,
# ||estimate - target||_F / ||target||_F:
#   theta_all, theta_unobs: the fitted matrix, on all entries and on the
#     unobserved ones;
#   beta_raw, gamma_raw: the estimates that the algorithm returns, against
#     the generating coefficients;
#   beta_id, gamma_id: the identified estimates, against their part of
#     T(theta).
# d_<name>: the relative difference between an estimate and that of the
# default-start fit of the same version, signal and threshold.
library(dplyr)
library(tidyr)
library(ggplot2)
library(RSSthemes)
options(width = 200)   # the tables are wide

size <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- if (identical(size, "lite")) "./output/sensitivity_lite" else
  "./output/sensitivity"
results <- readRDS(file.path(out_dir, "results_initialization.rds"))
cfg <- results$config

# Helpers ----
# "mean (standard deviation)" of a column, as text
mean_sd <- function(x, digits = 4) {
  return(sprintf("%.*f (%.*f)", digits, mean(x), digits, sd(x)))
}

# print a table under its title
show <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(as.data.frame(table), digits = 3, row.names = FALSE)
}

# the version of the truth and the kind of start, as labels
label_truth <- function(identified) {
  return(factor(identified, levels = c(TRUE, FALSE),
                labels = c("Truth satisfies (C1)–(C3)",
                           "Truth violates (C1)–(C3)")))
}
label_kind <- function(kind) {
  return(factor(kind, levels = c("default", "zero", "random"),
                labels = c("Default", "Zero", "Random")))
}

rows <- results$rows %>%
  mutate(truth = label_truth(identified), start = label_kind(kind)) %>%
  group_by(signal, identified, level) %>%
  mutate(
    # final objective, relative to the best start of the same version,
    # signal and threshold
    excess = (objective - min(objective)) / min(objective)
  ) %>%
  ungroup()

# Fits ----
show("Fits: iterations",
     rows %>% group_by(truth, level, start) %>%
       summarise(
         fits = n(), median_iterations = median(n_iter),
         largest = max(n_iter), not_converged = sum(!converged),
         rank = mean(rank_fit), seconds = mean(secs), .groups = "drop"
       ))

# Errors ----
errors <- c("theta_unobs", "beta_raw", "gamma_raw", "beta_id", "gamma_id")
show("Errors of the fits: mean (sd) over the fits",
     rows %>% group_by(truth, level, start) %>%
       summarise(across(all_of(errors), mean_sd), .groups = "drop"))

# Differences from the default start ----
# For each version of each signal: the largest difference from the
# default-start fit over the other starts, the largest excess of the final
# objective over the best start, and the range of the error on the
# unobserved entries over the starts.
differences <- c("d_theta", "d_beta_raw", "d_gamma_raw", "d_beta_id",
                 "d_gamma_id")
by_signal <- rows %>% group_by(level, truth, signal) %>%
  summarise(
    across(all_of(differences), max), excess = max(excess),
    theta_unobs_range = max(theta_unobs) - min(theta_unobs),
    .groups = "drop"
  )
show("Tight threshold, by signal: largest difference from the default start",
     by_signal %>% filter(level == "tight") %>% select(-level))

# The same, over the signals: the median and the largest of the differences
# of all the fits but the default-start ones.
show("Differences from the default start: median [largest] over the fits",
     rows %>% filter(kind != "default") %>% group_by(level, truth, start) %>%
       summarise(
         across(all_of(differences),
                ~ sprintf("%.1e [%.1e]", median(.x), max(.x))),
         .groups = "drop"
       ))

# Monotone objective ----
# An increase is a relative change of the objective above increase_tol from
# one iteration to the next.
show("Monotone objective",
     rows %>% group_by(truth) %>%
       summarise(
         fits = n(), iterations = sum(n_iter),
         with_an_increase = sum(n_increase),
         largest_relative_change = max(max_rel_increase), .groups = "drop"
       ))

# Figures ----
dir.create(file.path(out_dir, "paper"), showWarnings = FALSE)
save_figure <- function(figure, name, width, height) {
  ggsave(file.path(out_dir, "paper", paste0(name, ".pdf")), figure,
         device = cairo_pdf, width = width, height = height, units = "in")
}

# Figures of the paper ----
# Two figures for the first signal at the tight threshold, one line for each
# start: the objective after each iteration, and the errors at the
# checkpoints. As in identification_results.R, the labels use the symbols of
# the paper (the algorithm returns Theta hat, beta hat and Gamma hat), the
# theme and the colours are those of RSSthemes, and a figure is 8 in wide:
# reduced to the width of the text (about 6.4 in), its smallest text is
# about 8 pt.
theme_journal <- theme_significance() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(2.2, "lines"),
    legend.title = element_text(size = 12),
    legend.text = element_text(hjust = 0, size = 12),
    strip.text = element_text(face = "bold", size = 11),
    # the estimate, upright at the right of its row of panels
    strip.text.y = element_text(angle = 0, size = 13, face = "plain"),
    panel.spacing = unit(1.1, "lines")
  )
start_colors <- c(signif_blue, signif_orange, signif_red)

# The fits shown. The lines are drawn in the order of `fit`: the random
# starts first, then the zero start, and the default start on top.
shown <- function(data) {
  data <- data %>%
    filter(signal == cfg$signals[1], level == "tight") %>%
    mutate(truth = label_truth(identified), start = label_kind(kind)) %>%
    arrange(desc(start), start_id)
  return(data %>% mutate(fit = factor(fit, levels = unique(fit))))
}
path_lines <- list(
  geom_line(linewidth = 0.6, alpha = 0.85),
  scale_color_manual(values = start_colors),
  scale_x_log10(),
  labs(x = "Iteration (log scale)", color = "Initialization:"),
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1)))
)

# The objective after each iteration, minus the smallest value reached from
# any start in the same version of the signal.
traces <- shown(results$rows) %>%
  select(fit, truth, start) %>%
  left_join(
    bind_rows(lapply(names(results$objective), function(name) {
      objective <- results$objective[[name]]
      return(data.frame(fit = name, iteration = seq_along(objective),
                        objective = objective))
    })),
    by = "fit"
  ) %>%
  group_by(truth) %>%
  mutate(gap = objective - min(objective)) %>%
  ungroup() %>%
  filter(gap > 0)
objective_figure <- ggplot(traces, aes(iteration, gap, group = fit,
                                       color = start)) +
  path_lines +
  scale_y_log10(labels = scales::label_log()) +
  facet_wrap(~truth, nrow = 1) +
  labs(y = "Objective minus its smallest value\n(log scale)") +
  theme_journal
save_figure(objective_figure, "figure_sens_init_objective", width = 8,
            height = 4.2)

# The errors at the checkpoints: one row of panels for each estimate.
error_figure <- function(panels) {
  data <- shown(results$paths) %>%
    pivot_longer(all_of(names(panels)), names_to = "panel",
                 values_to = "error") %>%
    mutate(panel = factor(panel, levels = names(panels), labels = panels))
  figure <- ggplot(data, aes(iteration, error, group = fit, color = start)) +
    path_lines +
    scale_y_log10() +
    facet_grid(panel ~ truth, scales = "free_y",
               labeller = labeller(panel = label_parsed)) +
    labs(y = "RRMSE (log scale)") +
    theme_journal
  return(figure)
}
save_figure(
  error_figure(c(
    theta_unobs = 'atop(widehat(Theta) * ",", "unobserved")',
    beta_raw = "widehat(beta)", gamma_raw = "widehat(Gamma)"
  )),
  "figure_sens_init_rrmse", width = 8, height = 7.5
)

# Diagnostic figure (not in the paper) ----
# The identified estimates, against their part of T(Theta).
save_figure(
  error_figure(c(beta_id = "widetilde(beta)", gamma_id = "widetilde(Gamma)")),
  "figure_sens_init_rrmse_identified", width = 8, height = 5.6
)
