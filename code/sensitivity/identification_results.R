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
dir.create(file.path(out_dir, "paper"), showWarnings = FALSE)
save_figure <- function(figure, name, width, height) {
  ggsave(file.path(out_dir, "paper", paste0(name, ".pdf")), figure,
         device = cairo_pdf, width = width, height = height, units = "in")
}

# Figures of the paper ----
# Three figures: the errors against the SNR, the errors against the missing
# rate, and the followed coefficients of one signal. Their labels use the
# symbols of the paper: the estimate that the algorithm returns is beta hat
# (Gamma hat), the identified one is beta tilde (Gamma tilde), a true value
# has the superscript (*), and T(Theta^(*)) holds the identified version of
# the true coefficients. In every figure of the paper, beta hat has a solid
# line and filled circles, and beta tilde a dashed line and open triangles.
# The theme and the colours are those of RSSthemes, as in the other figures
# of the paper. A figure is 8 in wide: reduced to the width of the text
# (about 6.4 in), its smallest text is about 8 pt.
theme_journal <- theme_significance() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.box = "vertical", legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(1, "pt"),
    legend.key.width = unit(2.6, "lines"),
    legend.text = element_text(hjust = 0, size = 12),
    strip.text = element_text(face = "bold", size = 11),
    # beta and Gamma, upright at the right of their row of panels
    strip.text.y = element_text(angle = 0, size = 15),
    panel.spacing = unit(1.1, "lines")
  )
estimate_labels <- expression(
  widehat(beta) * "," ~ widehat(Gamma),
  widetilde(beta) * "," ~ widetilde(Gamma)
)
estimate_lines <- c("solid", "22")
estimate_shapes <- c(16, 2)

# Error on beta and on Gamma against the level of a scenario: mean and
# +/- 1 standard deviation over the signals. The line type and the shape
# give the estimate, the colour gives what it is compared with: the true
# coefficients or their identified version. With a truth that satisfies the
# conditions the two are the same, and only the first is drawn.
target_labels <- expression(
  "the true coefficients" ~ beta^"(*)" * "," ~ Gamma^"(*)",
  "their identified version in" ~ T(Theta^"(*)")
)
# an orange and a dark blue: the two remain distinct when printed in grey
target_colors <- c(signif_orange, "#00485b")
curves <- by_signal %>%
  pivot_longer(matches("^(beta|gamma)_(id|raw)_(T|gen)$"),
               names_to = c("component", "estimate", "target"),
               names_pattern = "^(beta|gamma)_(id|raw)_(T|gen)$",
               values_to = "error") %>%
  filter(!(identified & target == "T")) %>%
  group_by(scenario, truth, component, estimate, target, level) %>%
  summarise(mean = mean(error), sd = sd(error), .groups = "drop") %>%
  mutate(estimate = factor(estimate, levels = c("raw", "id")),
         target = factor(target, levels = c("gen", "T")),
         component = factor(component, levels = c("beta", "gamma"),
                            labels = c("beta", "Gamma")))

error_figure <- function(scenario_name, x_label, x_scale) {
  data <- curves %>% filter(scenario == scenario_name)
  figure <- ggplot(data, aes(level, mean, color = target, fill = target,
                             linetype = estimate, shape = estimate,
                             group = interaction(estimate, target))) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.12,
                color = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.3, stroke = 0.9) +
    scale_color_manual(values = target_colors, labels = target_labels) +
    scale_fill_manual(values = target_colors, guide = "none") +
    scale_linetype_manual(values = estimate_lines, labels = estimate_labels) +
    scale_shape_manual(values = estimate_shapes, labels = estimate_labels) +
    # two legends, one under the other: the estimate, then the comparison
    guides(
      linetype = guide_legend(order = 1), shape = guide_legend(order = 1),
      color = guide_legend(order = 2, override.aes = list(shape = NA))
    ) +
    x_scale(breaks = sort(unique(data$level))) +
    expand_limits(y = 0) +
    facet_grid(component ~ truth,
               labeller = labeller(component = label_parsed)) +
    labs(x = x_label, y = "RRMSE", linetype = "Estimate:",
         shape = "Estimate:", color = "RRMSE against:") +
    theme_journal
  return(figure)
}
save_figure(error_figure("Noise", "Signal-to-noise ratio (log scale)",
                         scale_x_log10),
            "figure_sens_ident_noise", width = 8, height = 6)
save_figure(error_figure("Missingness", "Missing entries (%)",
                         scale_x_continuous),
            "figure_sens_ident_missing", width = 8, height = 6)

# The followed coefficients of the first signal in the cell (snr_base,
# missing_base). For each coefficient: the two estimates minus the true
# value, in the units of the coefficients, with +/- 1 standard deviation
# over the replications. The cross is where the identified version of the
# true coefficient is: at zero when the truth satisfies the conditions.
# The coefficients are ordered by the height of that cross in the truth
# that violates the conditions, in both columns of panels.
base_coefs <- by_coef %>%
  filter(signal == cfg$signals[1], snr == cfg$snr_base,
         missing == cfg$missing_base)
coef_order <- base_coefs %>%
  filter(!identified, estimate == "Raw") %>%
  group_by(component) %>%
  mutate(position = rank(truth_T, ties.method = "first")) %>%
  ungroup() %>%
  select(component, coef, position)
# one row for each mark: the two estimates (with their standard deviation)
# and the cross
mark_labels <- c(
  estimate_labels,
  expression("identified version of the true coefficient, in" ~
               T(Theta^"(*)"))
)
marks <- bind_rows(
  base_coefs %>%
    mutate(mark = if_else(estimate == "Raw", "raw", "id"), value = error_gen),
  base_coefs %>%
    filter(estimate == "Raw") %>%
    mutate(mark = "target", value = truth_T, sd = NA_real_)
) %>%
  left_join(coef_order, by = c("component", "coef")) %>%
  mutate(
    # the estimates on each side of the cross
    x = position + c(raw = -0.22, target = 0, id = 0.22)[mark],
    mark = factor(mark, levels = c("raw", "id", "target")),
    truth = label_truth(identified),
    component = factor(component, levels = c("beta", "gamma"),
                       labels = c("beta", "Gamma"))
  )
coef_figure <- ggplot(marks, aes(x, value, color = mark, shape = mark)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_linerange(data = marks %>% filter(mark != "target"),
                 aes(ymin = value - sd, ymax = value + sd), linewidth = 0.4,
                 show.legend = FALSE) +
  geom_point(size = 1.9, stroke = 0.8) +
  scale_color_manual(values = c(signif_red, signif_blue, "black"),
                     labels = mark_labels) +
  scale_shape_manual(values = c(estimate_shapes, 4), labels = mark_labels) +
  scale_x_continuous(breaks = seq_len(max(marks$position)), labels = NULL) +
  facet_grid(component ~ truth,
             labeller = labeller(component = label_parsed)) +
  labs(x = "Coefficient", y = "Difference from the true coefficient",
       color = NULL, shape = NULL) +
  theme_journal +
  theme(panel.grid.minor.x = element_blank())
save_figure(coef_figure, "figure_sens_ident_coefs", width = 8, height = 5.6)

# Diagnostic figure (not in the paper) ----
theme_paper <- theme_significance() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold", size = 11))

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
