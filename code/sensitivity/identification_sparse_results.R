# Sensitivity analysis, identified components with sparse coefficients:
# tables and figures.
#
# Reads results_identification_sparse.rds (written by
# identification_sparse_run.R), prints the tables and saves the figures in
# the folder paper/ of the run. It follows identification_results.R.
# Run from the root of the simulation repository, with `lite` as an argument
# for the test run:
#   Rscript code/sensitivity/identification_sparse_results.R
#
# One row of the results is one fit, with lambda_beta, lambda_gamma and
# lambda_m tuned. An error is the RRMSE of the paper,
# ||estimate - target||_F / ||target||_F. The errors are named
# <component>_<estimate>_<target>:
#   component: beta, gamma or M;
#   estimate : raw (as the algorithm returns it), id (identified) or plain
#              (raw, from the fit without Lasso penalties on the same data);
#   target   : gen (the generating component, which is sparse) or T (its
#              part of T(theta), which is not).
# The truth is the one as generated: it violates (C1)-(C3) by chance only.
#
# The fits of a signal are first averaged over its replications. The tables
# and the figures then give the mean over the signals, with its standard
# error or the standard deviation over the signals.
library(dplyr)
library(tidyr)
library(ggplot2)
library(RSSthemes)

size <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- if (identical(size, "lite")) "./output/sensitivity_lite" else
  "./output/sensitivity"
results <- readRDS(file.path(out_dir, "results_identification_sparse.rds"))
cfg <- results$config

# Helpers ----
# "mean (standard error)" of a column, as text
mean_se <- function(x, digits = 4) {
  return(sprintf("%.*f (%.*f)", digits, mean(x), digits,
                 sd(x) / sqrt(length(x))))
}

# "mean [smallest, largest]" of a column, as text
mean_range <- function(x, digits = 3) {
  return(sprintf("%.*f [%.*f, %.*f]", digits, mean(x), digits, min(x),
                 digits, max(x)))
}

# print a table under its title
show <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(as.data.frame(table), digits = 3, row.names = FALSE)
}

# The two scenarios, one under the other. `level` is the SNR in the scenario
# "Noise" and the percentage of missing entries in the scenario
# "Missingness". The cell (snr_base, missing_base) is in both.
stack_scenarios <- function(data) {
  stacked <- bind_rows(
    data %>% filter(missing == cfg$missing_base) %>%
      mutate(scenario = "Noise", level = snr),
    data %>% filter(snr == cfg$snr_base) %>%
      mutate(scenario = "Missingness", level = 100 * missing)
  )
  return(stacked %>% mutate(
    scenario = factor(scenario, levels = c("Noise", "Missingness"))
  ))
}

# Signals ----
# gap: the RRMSE of a generating component against its part of T(theta).
# The identified estimate cannot be closer than that to the generating one.
show("The signals: mean [smallest, largest] over the signals",
     results$truth$facts %>%
       summarise(
         signals = n(),
         across(c(sd_theta, gap_beta, gap_gamma, gap_M), mean_range)
       ))

# Fits ----
# One row per signal and cell: the means over the replications of the
# signal. The tables below are means over the signals.
by_signal <- stack_scenarios(results$rows) %>%
  group_by(scenario, level, signal) %>%
  summarise(across(c(where(is.numeric), converged), mean), fits = n(),
            .groups = "drop")
by_cell <- by_signal %>% group_by(scenario, level)

# lambda_*_share: the selected penalty over the top of its grid. A share of
# 0 or 1 means that the grid should be moved.
show("Fits and selected penalties",
     by_cell %>% summarise(
       signals = n(), fits = sum(fits), converged = mean(converged),
       lambda_beta = mean(lambda_beta),
       lambda_beta_share = mean(lambda_beta / lambda_beta_max),
       lambda_gamma = mean(lambda_gamma),
       lambda_gamma_share = mean(lambda_gamma / lambda_gamma_max),
       lambda_m = mean(lambda_m), sigma_hat = mean(sigma_hat),
       sigma = mean(sigma), seconds = mean(secs), .groups = "drop"
     ))

show("Selected penalties at an end of their grid: share of the fits",
     stack_scenarios(results$rows) %>% group_by(scenario, level) %>%
       summarise(
         beta_at_0 = mean(lambda_beta == 0),
         beta_at_top = mean(lambda_beta > 0.999 * lambda_beta_max),
         gamma_at_0 = mean(lambda_gamma == 0),
         gamma_at_top = mean(lambda_gamma > 0.999 * lambda_gamma_max),
         .groups = "drop"
       ))

# Errors ----
# theta_all_plain: the fit without Lasso penalties.
show("Fitted matrix (with Lasso, without, unobserved entries, 1 - R^2)",
     by_cell %>% summarise(
       across(c(theta_all, theta_all_plain, theta_unobs, one_minus_r2),
              mean_se),
       .groups = "drop"
     ))

show("beta: mean (SE)",
     by_cell %>% summarise(
       across(c(beta_raw_gen, beta_plain_gen, beta_id_gen, beta_raw_T,
                beta_id_T), mean_se),
       .groups = "drop"
     ))

show("Gamma: mean (SE)",
     by_cell %>% summarise(
       across(c(gamma_raw_gen, gamma_plain_gen, gamma_id_gen, gamma_raw_T,
                gamma_id_T), mean_se),
       .groups = "drop"
     ))

show("M: mean (SE)",
     by_cell %>% summarise(
       across(c(M_raw_gen, M_id_gen, M_raw_T, M_id_T), mean_se),
       .groups = "drop"
     ))

# Zeros ----
# *_zeros: share of the entries of the raw estimate that are exactly zero
# (half of the true ones are). *_zeros_found: share of the true zeros
# estimated as zero. *_nonzeros_kept: share of the true non-zeros estimated
# as non-zero.
show("Zeros of the raw estimates: means",
     by_cell %>% summarise(
       across(c(beta_zeros, beta_zeros_found, beta_nonzeros_kept,
                gamma_zeros, gamma_zeros_found, gamma_nonzeros_kept), mean),
       .groups = "drop"
     ))

# Raw against identified ----
# *_raw_id: the distance between the two estimates. c1, c2, c3: how far the
# raw split is from (C1) X'M = 0, (C2) MZ = 0 and (C3) X'Gamma = 0.
show("Raw against identified: means",
     by_cell %>% summarise(
       beta_raw_id = mean(beta_raw_id), gamma_raw_id = mean(gamma_raw_id),
       M_raw_id = mean(M_raw_id),
       c1 = mean(c1_raw), c2 = mean(c2_raw), c3 = mean(c3_raw),
       .groups = "drop"
     ))

# Each of these is zero up to rounding error when the code is right: the
# change of the fitted matrix by the identification step, the sum of the
# three shares of its error minus one, and the distance of the identified
# split from (C1)-(C3).
show("Checks: largest value over all fits",
     results$rows %>% summarise(
       theta_changed = max(theta_same),
       shares_minus_one = max(abs(share_beta + share_gamma + share_M - 1)),
       c1_id = max(c1_id), c2_id = max(c2_id), c3_id = max(c3_id)
     ))

# Followed coefficients ----
# As in identification_results.R: the mean and the standard deviation of
# each followed estimate over the replications, in % of the typical size of
# a coefficient (the root mean square of the entries of the true matrix).
# `zero`: the generating value is zero.
by_coef <- results$coefs %>%
  left_join(results$truth$coefs, by = c("signal", "identified", "coef")) %>%
  pivot_longer(c(raw, id), names_to = "estimate", values_to = "value") %>%
  group_by(signal, snr, missing, component, coef, estimate, gen, T_theta,
           size_gen, size_T) %>%
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(
    error_gen = 100 * (mean - gen) / size_gen,
    error_T = 100 * (mean - T_theta) / size_T,
    sd = 100 * sd / size_gen,
    truth_T = 100 * (T_theta - gen) / size_gen,
    zero = gen == 0,
    estimate = factor(estimate, levels = c("raw", "id"),
                      labels = c("Raw", "Identified"))
  )

# Means over the followed coefficients of all the signals, for the true
# zeros and for the others.
show("Followed coefficients, in % of the typical size of a coefficient",
     stack_scenarios(by_coef) %>%
       group_by(scenario, level, component, zero, estimate) %>%
       summarise(
         abs_error_gen = mean(abs(error_gen)),
         abs_error_T = mean(abs(error_T)), sd = mean(sd),
         coefficients = n(), .groups = "drop"
       ))

# Figures ----
# The theme and the colours are those of identification_results.R.
theme_paper <- theme_significance() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold", size = 11))
legend_rows <- guide_legend(nrow = 2, byrow = TRUE)
dir.create(file.path(out_dir, "paper"), showWarnings = FALSE)
save_figure <- function(figure, name, width, height) {
  ggsave(file.path(out_dir, "paper", paste0(name, ".pdf")), figure,
         device = cairo_pdf, width = width, height = height, units = "in")
}
component_labels <- c(beta = "beta", gamma = "Gamma")

# Error on beta and on Gamma against the level of a scenario: mean and
# +/- 1 standard deviation over the signals.
line_labels <- c(
  id_T = "Identified estimate vs. identified truth",
  id_gen = "Identified estimate vs. generating coefficients",
  raw_gen = "Raw estimate vs. generating coefficients",
  plain_gen = "Raw estimate without Lasso vs. generating coefficients"
)
line_colors <- setNames(
  c(signif_blue, signif_orange, signif_red, "grey40"), line_labels
)
line_shapes <- setNames(c(16, 17, 15, 4), line_labels)
curves <- by_signal %>%
  pivot_longer(matches("^(beta|gamma)_(id_T|id_gen|raw_gen|plain_gen)$"),
               names_to = c("component", "line"),
               names_pattern = "^(beta|gamma)_(.*)$", values_to = "error") %>%
  group_by(scenario, component, line, level) %>%
  summarise(mean = mean(error), sd = sd(error), .groups = "drop") %>%
  mutate(line = factor(line, levels = names(line_labels),
                       labels = line_labels),
         component = factor(component, levels = names(component_labels),
                            labels = component_labels))

error_figure <- function(scenario_name, x_label, x_scale) {
  data <- curves %>% filter(scenario == scenario_name)
  figure <- ggplot(data, aes(level, mean, color = line, fill = line,
                             shape = line)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15,
                color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = line_colors) +
    scale_fill_manual(values = line_colors) +
    scale_shape_manual(values = line_shapes) +
    # the legend in two rows: the identified estimate, then the raw one
    guides(color = legend_rows, fill = legend_rows, shape = legend_rows) +
    x_scale(breaks = sort(unique(data$level))) +
    facet_wrap(~component, nrow = 1, labeller = label_parsed) +
    labs(x = x_label, y = "RRMSE", color = NULL, fill = NULL, shape = NULL) +
    theme_paper
  return(figure)
}
save_figure(error_figure("Noise", "Signal-to-noise ratio (log scale)",
                         scale_x_log10),
            "figure_sens_ident_sparse_noise", width = 12, height = 5)
save_figure(error_figure("Missingness", "Missing entries (%)",
                         scale_x_continuous),
            "figure_sens_ident_sparse_missing", width = 12, height = 5)

# The zeros of the raw estimates against the level of a scenario: the share
# of the true zeros estimated as zero, and the share of the true non-zeros
# estimated as non-zero. Means over the signals.
zero_labels <- c(zeros_found = "True zeros estimated as zero",
                 nonzeros_kept = "True non-zeros estimated as non-zero")
zeros <- by_signal %>%
  pivot_longer(matches("^(beta|gamma)_(zeros_found|nonzeros_kept)$"),
               names_to = c("component", "line"),
               names_pattern = "^(beta|gamma)_(.*)$", values_to = "share") %>%
  group_by(scenario, component, line, level) %>%
  summarise(mean = mean(share), .groups = "drop") %>%
  mutate(line = factor(line, levels = names(zero_labels),
                       labels = zero_labels),
         component = factor(component, levels = names(component_labels),
                            labels = component_labels))

zeros_figure <- function(scenario_name, x_label, x_scale) {
  data <- zeros %>% filter(scenario == scenario_name)
  figure <- ggplot(data, aes(level, mean, color = line, shape = line)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c(signif_blue, signif_orange)) +
    scale_shape_manual(values = c(16, 17)) +
    x_scale(breaks = sort(unique(data$level))) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~component, nrow = 1, labeller = label_parsed) +
    labs(x = x_label, y = "Share", color = NULL, shape = NULL) +
    theme_paper
  return(figure)
}
save_figure(zeros_figure("Noise", "Signal-to-noise ratio (log scale)",
                         scale_x_log10),
            "figure_sens_ident_sparse_zeros_noise", width = 12, height = 5)
save_figure(zeros_figure("Missingness", "Missing entries (%)",
                         scale_x_continuous),
            "figure_sens_ident_sparse_zeros_missing", width = 12, height = 5)

# The followed coefficients of the first signal in the cell (snr_base,
# missing_base). For each coefficient: the error against the generating
# value, in % of the typical size of a coefficient, with +/- 1 standard
# deviation over the replications. The cross is where T(theta) is. The
# left panels hold the true zeros, the right ones the other coefficients.
base_coefs <- by_coef %>%
  filter(signal == cfg$signals[1], snr == cfg$snr_base,
         missing == cfg$missing_base) %>%
  group_by(component, zero) %>%
  mutate(position = match(coef, sort(unique(coef)))) %>%
  ungroup() %>%
  mutate(kind = factor(zero, levels = c(TRUE, FALSE),
                       labels = c("True value is zero",
                                  "True value is not zero")),
         component = factor(component, levels = names(component_labels),
                            labels = component_labels))
dodge <- position_dodge(width = 0.6)
coef_figure <- ggplot(base_coefs, aes(position, error_gen, color = estimate,
                                      shape = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(aes(y = truth_T), color = "black", shape = 4, size = 2) +
  geom_errorbar(aes(ymin = error_gen - sd, ymax = error_gen + sd), width = 0,
                position = dodge) +
  geom_point(size = 2, position = dodge) +
  scale_color_manual(values = c(signif_red, signif_blue)) +
  scale_shape_manual(values = c(15, 16)) +
  facet_grid(component ~ kind, scales = "free",
             labeller = labeller(component = label_parsed)) +
  labs(x = "Followed coefficient",
       y = "Estimate minus true value (% of the typical coefficient size)",
       color = NULL, shape = NULL) +
  theme_paper
save_figure(coef_figure, "figure_sens_ident_sparse_coefs", width = 12,
            height = 7)
