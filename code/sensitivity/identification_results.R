# Sensitivity analysis, identified components: tables and figures.
#
# Reads results_identification.rds (written by identification_run.R), prints
# the tables and saves the figures in the folder paper/ of the run.
# Run from the root of the simulation repository:
#   Rscript code/sensitivity/identification_results.R         (draft or full)
#   Rscript code/sensitivity/identification_results.R lite    (the test run)
#
# One row of the results is one fit. An error is the RRMSE of the paper,
# ||estimate - target||_F / ||target||_F. The errors are named
# <component>_<estimate>_<target>:
#   component: beta, gamma or M;
#   estimate : raw (as the algorithm returns it) or id (identified);
#   target   : gen (the generating component) or T (its part of T(theta)).
# <component>_raw_id is the distance between the two estimates.
# one_minus_r2 is 1 - R^2 of the fitted matrix: its squared error over the
# variation of theta around its mean.
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
results <- readRDS(file.path(out_dir, "results_identification.rds"))
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

# the version of the truth, as a label
label_truth <- function(identified) {
  return(factor(identified, levels = c(TRUE, FALSE),
                labels = c("Truth satisfies (C1)–(C3)",
                           "Truth violates (C1)–(C3)")))
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
    scenario = factor(scenario, levels = c("Noise", "Missingness")),
    truth = label_truth(identified)
  ))
}

# Signals ----
# gap: the RRMSE of a generating component against its part of T(theta).
# The identified estimate cannot be closer than that to the generating one.
show("The signals: mean [smallest, largest] over the signals",
     results$truth$facts %>% group_by(identified) %>%
       summarise(
         signals = n(),
         across(c(sd_theta, gap_beta, gap_gamma, gap_M), mean_range),
         .groups = "drop"
       ))

# Fits ----
# One row per signal and cell: the means over the replications of the
# signal. The tables below are means over the signals.
by_signal <- stack_scenarios(results$rows) %>%
  group_by(scenario, truth, identified, level, signal) %>%
  summarise(across(c(where(is.numeric), converged), mean), fits = n(),
            .groups = "drop")
by_cell <- by_signal %>% group_by(scenario, truth, level)

show("Fits",
     by_cell %>% summarise(
       signals = n(), fits = sum(fits), converged = mean(converged),
       iterations = mean(n_iter), lambda_m = mean(lambda_m),
       rank = mean(rank_fit), seconds = mean(secs), .groups = "drop"
     ))

# Errors ----
show("Fitted matrix (all entries, unobserved entries, 1 - R^2): mean (SE)",
     by_cell %>% summarise(
       across(c(theta_all, theta_unobs, one_minus_r2), mean_se),
       .groups = "drop"
     ))

show("M: mean (SE)",
     by_cell %>% summarise(
       across(c(M_raw_gen, M_id_gen, M_raw_T, M_id_T), mean_se),
       .groups = "drop"
     ))

show("beta: mean (SE)",
     by_cell %>% summarise(
       across(c(beta_raw_gen, beta_id_gen, beta_raw_T, beta_id_T), mean_se),
       .groups = "drop"
     ))

show("Gamma: mean (SE)",
     by_cell %>% summarise(
       across(c(gamma_raw_gen, gamma_id_gen, gamma_raw_T, gamma_id_T),
              mean_se),
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

# The squared error of the fitted matrix is the sum of three parts, one for
# each identified component.
show("Shares of the squared error of the fitted matrix: means",
     by_cell %>% summarise(
       beta = mean(share_beta), gamma = mean(share_gamma), M = mean(share_M),
       .groups = "drop"
     ))

# Each of these is zero up to rounding error when the code is right: the
# change of the fitted matrix by the identification step, the sum of the
# three shares minus one, and the distance of the identified split from
# (C1)-(C3).
show("Checks: largest value over all fits",
     results$rows %>% summarise(
       theta_changed = max(theta_same),
       shares_minus_one = max(abs(share_beta + share_gamma + share_M - 1)),
       c1_id = max(c1_id), c2_id = max(c2_id), c3_id = max(c3_id)
     ))

# Followed coefficients ----
# For each followed coefficient of each signal, in each cell and for each
# estimate: the mean and the standard deviation of the estimate over the
# replications, in the units of the coefficients (nothing is scaled).
#   error_gen: mean estimate - generating value;
#   error_T  : mean estimate - value in T(theta);
#   sd       : sd(estimate);
#   truth_T  : where T(theta) is, T - gen.
by_coef <- results$coefs %>%
  left_join(results$truth$coefs, by = c("signal", "identified", "coef")) %>%
  pivot_longer(c(raw, id), names_to = "estimate", values_to = "value") %>%
  group_by(signal, identified, snr, missing, component, coef, estimate, gen,
           T_theta) %>%
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(
    error_gen = mean - gen,
    error_T = mean - T_theta,
    truth_T = T_theta - gen,
    estimate = factor(estimate, levels = c("raw", "id"),
                      labels = c("Raw", "Identified"))
  )

# For the scale of the errors below: the size of the true coefficients.
show("Followed coefficients: size of the true values",
     results$truth$coefs %>%
       mutate(truth = label_truth(identified)) %>%
       group_by(truth, component) %>%
       summarise(
         coefficients = n(), mean_abs = mean(abs(gen)),
         root_mean_square = sqrt(mean(gen^2)), .groups = "drop"
       ))

# Means over the followed coefficients of all the signals.
show("Followed coefficients: errors in the units of the coefficients",
     stack_scenarios(by_coef) %>%
       group_by(scenario, truth, level, component, estimate) %>%
       summarise(
         abs_error_gen = mean(abs(error_gen)),
         abs_error_T = mean(abs(error_T)), sd = mean(sd),
         coefficients = n(), .groups = "drop"
       ))

# Figures ----
# The theme and the colours are those of 4_plots.R.
theme_paper <- theme_significance() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold", size = 11))
legend_rows <- guide_legend(nrow = 2, byrow = TRUE)
dir.create(file.path(out_dir, "paper"), showWarnings = FALSE)
save_figure <- function(figure, name, width, height) {
  ggsave(file.path(out_dir, "paper", paste0(name, ".pdf")), figure,
         device = cairo_pdf, width = width, height = height, units = "in")
}

# Error on beta and on Gamma against the level of a scenario: mean and
# +/- 1 standard deviation over the signals. One line for each estimate
# (identified, raw) against each target (identified truth, generating
# coefficients). With a truth that satisfies the conditions the two targets
# are the same, so each estimate has one line there.
line_labels <- c(
  id_T = "Identified estimate vs. identified truth",
  id_gen = "Identified estimate vs. generating coefficients",
  raw_T = "Raw estimate vs. identified truth",
  raw_gen = "Raw estimate vs. generating coefficients"
)
line_colors <- setNames(
  c(signif_blue, signif_orange, signif_green, signif_red), line_labels
)
line_shapes <- setNames(c(16, 17, 18, 15), line_labels)
curves <- by_signal %>%
  pivot_longer(matches("^(beta|gamma)_(id|raw)_(T|gen)$"),
               names_to = c("component", "line"),
               names_pattern = "^(beta|gamma)_(.*)$", values_to = "error") %>%
  filter(!(identified & line %in% c("id_gen", "raw_T"))) %>%
  group_by(scenario, truth, component, line, level) %>%
  summarise(mean = mean(error), sd = sd(error), .groups = "drop") %>%
  mutate(line = factor(line, levels = names(line_labels),
                       labels = line_labels),
         component = factor(component, levels = c("beta", "gamma"),
                            labels = c("beta", "Gamma")))

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
    facet_grid(component ~ truth, scales = "free_y",
               labeller = labeller(component = label_parsed)) +
    labs(x = x_label, y = "RRMSE", color = NULL, fill = NULL, shape = NULL) +
    theme_paper
  return(figure)
}
save_figure(error_figure("Noise", "Signal-to-noise ratio (log scale)",
                         scale_x_log10),
            "figure_sens_ident_noise", width = 12, height = 7)
save_figure(error_figure("Missingness", "Missing entries (%)",
                         scale_x_continuous),
            "figure_sens_ident_missing", width = 12, height = 7)

# The followed coefficients of the first signal in the cell (snr_base,
# missing_base). For each coefficient: the error against the generating
# value, in the units of the coefficients, with +/- 1 standard deviation
# over the replications. The cross is where T(theta) is: at zero when the
# truth satisfies the conditions.
base_coefs <- by_coef %>%
  filter(signal == cfg$signals[1], snr == cfg$snr_base,
         missing == cfg$missing_base) %>%
  group_by(component) %>%
  mutate(position = coef - min(coef) + 1) %>%
  ungroup() %>%
  mutate(truth = label_truth(identified),
         component = factor(component, levels = c("beta", "gamma"),
                            labels = c("beta", "Gamma")))
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
  facet_grid(component ~ truth, scales = "free_y",
             labeller = labeller(component = label_parsed)) +
  labs(x = "Followed coefficient", y = "Estimate minus true value",
       color = NULL, shape = NULL) +
  theme_paper
save_figure(coef_figure, "figure_sens_ident_coefs", width = 12, height = 7)

# Every followed coefficient of every signal in the cell (snr_base,
# missing_base): the error of its estimate, averaged over the replications,
# against its true value. No coefficient is left out and nothing is scaled.
all_coefs <- by_coef %>%
  filter(snr == cfg$snr_base, missing == cfg$missing_base) %>%
  mutate(truth = label_truth(identified),
         component = factor(component, levels = c("beta", "gamma"),
                            labels = c("beta", "Gamma")))
all_figure <- ggplot(all_coefs, aes(gen, error_gen, color = estimate,
                                    shape = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_manual(values = c(signif_red, signif_blue)) +
  scale_shape_manual(values = c(15, 16)) +
  facet_grid(component ~ truth + estimate,
             labeller = labeller(component = label_parsed)) +
  labs(x = "True value", y = "Estimate minus true value") +
  theme_paper +
  theme(legend.position = "none")
save_figure(all_figure, "figure_sens_ident_coefs_all", width = 12,
            height = 7)
