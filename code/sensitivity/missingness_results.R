# Sensitivity analysis, missingness mechanism: tables and figure.
#
# Reads results_missingness.rds (written by missingness_run.R), prints the
# tables and saves the figure in the folder paper/ of the run.
# Run from the root of the simulation repository:
#   Rscript code/sensitivity/missingness_results.R         (draft or full run)
#   Rscript code/sensitivity/missingness_results.R lite    (the test run)
#
# One row of the results is one fit. All errors are relative,
# ||estimate - truth||_F / ||truth||_F. The first rows are the rows that
# every size contains.
#   theta_all, theta_unobs: the fitted matrix on the first rows, on all
#     their entries and on the unobserved ones;
#   theta_low, theta_high: the same on all entries of the tenth of those
#     rows that the MAR mechanisms observe least, and of the tenth that they
#     observe most;
#   theta_unobs_full: the fitted matrix on the unobserved entries of all the
#     rows of the fit;
#   beta_raw_t, beta_id_t: the raw and the identified estimates of beta,
#     against the identified truth of the rows of the fit, T(theta). It is
#     the generated beta at the largest size only;
#   beta_raw, beta_id: the same estimates against the generated beta;
#   gamma_raw_t, gamma_id_t, gamma_raw, gamma_id: the same four errors for
#     Gamma, on the first rows;
#   gamma_raw_full, gamma_id_full: the raw and the identified Gamma against
#     the generated Gamma, on all the rows of the fit;
#   gamma_low, gamma_high: the raw Gamma against the identified truth, on
#     the least and on the most observed tenth of the first rows.
library(dplyr)
library(tidyr)
library(ggplot2)
library(RSSthemes)
options(width = 200)   # the tables are wide

size <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- if (identical(size, "lite")) "./output/sensitivity_lite" else
  "./output/sensitivity"
rows <- readRDS(file.path(out_dir, "results_missingness.rds"))$rows

# Helpers ----
# "mean (Monte Carlo standard error)" of a column, as text
mean_se <- function(x, digits = 4) {
  sprintf("%.*f (%.*f)", digits, mean(x), digits, sd(x) / sqrt(length(x)))
}

# print a table under its title
show <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(as.data.frame(table), digits = 3, row.names = FALSE)
}

# Fits ----
# A fit is skipped when its mask leaves a row or column with fewer than two
# observed entries.
show("Fits by mechanism and status",
     rows %>% count(kappa_x, kappa_z, status, name = "fits"))

# The comparisons are paired: only the replications whose fits all ran.
rows <- rows %>% group_by(seed) %>% filter(all(status == "ok")) %>% ungroup()

# Mechanisms ----
# logit(prob_ij) = a0 + kappa_x * score_i + kappa_z * score_j. Both zero:
# MCAR. "MAR X": driven by X only (kappa_z = 0). "MAR XZ": driven by X and Z
# equally (kappa_z = kappa_x).
# pi_ratio: the largest probability of observation over the smallest, and
# obs_rate: the observed rate, both on the largest size. The fewest entries
# in a row and in a column are taken over all the sizes.
mechanisms <- rows %>% group_by(kappa_x, kappa_z) %>%
  summarise(
    pi_ratio = mean(pi_ratio[n == max(n)]),
    obs_rate = mean(obs_rate[n == max(n)]),
    fewest_in_a_row = min(min_row_obs),
    fewest_in_a_column = min(min_col_obs),
    replications = n_distinct(seed), .groups = "drop"
  ) %>%
  arrange(kappa_z > 0, kappa_x) %>%
  mutate(mechanism = case_when(
    kappa_x == 0 ~ "MCAR",
    kappa_z == 0 ~ sprintf("MAR X (%.0f)", pi_ratio),
    TRUE ~ sprintf("MAR XZ (%.0f)", pi_ratio)
  ))
show("Mechanisms", mechanisms)

rows <- rows %>%
  left_join(select(mechanisms, kappa_x, kappa_z, mechanism),
            by = c("kappa_x", "kappa_z")) %>%
  mutate(mechanism = factor(mechanism, levels = mechanisms$mechanism))

# Errors ----
# Two tables for each group of errors: the errors, then the errors minus
# their MCAR value in the same replication and at the same n. The SE of a
# difference is small because the fits are paired.
error_groups <- list(
  "Fitted matrix" = c("theta_all", "theta_unobs", "theta_unobs_full",
                      "theta_low", "theta_high"),
  "Beta" = c("beta_raw_t", "beta_id_t", "beta_raw", "beta_id"),
  "Gamma" = c("gamma_raw_t", "gamma_id_t", "gamma_raw", "gamma_id",
              "gamma_raw_full", "gamma_id_full", "gamma_low", "gamma_high")
)
for (group in names(error_groups)) {
  errors <- error_groups[[group]]
  show(paste0(group, ", errors: mean (SE)"),
       rows %>% group_by(mechanism, n) %>%
         summarise(across(all_of(errors), mean_se), .groups = "drop"))
  show(paste0(group, ", difference from MCAR within a replication: mean (SE)"),
       rows %>% group_by(seed, n) %>%
         mutate(across(all_of(errors), ~ .x - .x[kappa_x == 0])) %>%
         ungroup() %>% filter(kappa_x > 0) %>% group_by(mechanism, n) %>%
         summarise(across(all_of(errors), mean_se), .groups = "drop"))
}

# Tuning and checks ----
# beta_raw_id, gamma_raw_id: relative distance between the raw and the
#   identified estimates (Gamma: first rows).
# beta_gap, gamma_gap: relative distance between the generated beta and
#   Gamma and the identified truth of the first n rows, T(theta). They are
#   zero at the largest size, where conditions (C1)-(C3) hold. At the other
#   sizes the estimates target that truth, so these distances are part of
#   their errors against the generated coefficients.
show("Tuning and checks: means",
     rows %>% group_by(mechanism, n) %>%
       summarise(
         lambda_m = mean(lambda_m), rank = mean(rank_fit),
         iterations = mean(n_iter), converged = mean(converged),
         seconds = mean(secs),
         beta_raw_id = mean(beta_raw_id), gamma_raw_id = mean(gamma_raw_id),
         beta_gap = mean(beta_gap), gamma_gap = mean(gamma_gap),
         .groups = "drop"
       ))

# Figure ----
# Error against n, one line per mechanism, for every error of the tables.
# The panels of the paper are chosen from these later. The theme is that of
# 4_plots.R. MCAR is grey; the mechanisms driven by X are three steps of
# blue (triangles) and those driven by X and Z three steps of orange
# (squares), from the weakest to the strongest.
panels <- c(
  theta_all = "Fitted matrix: all entries",
  theta_unobs = "Fitted matrix: unobserved entries",
  theta_low = "Fitted matrix: least observed rows",
  theta_high = "Fitted matrix: most observed rows",
  beta_raw_t = "Raw beta vs identified truth",
  beta_id_t = "Identified beta vs identified truth",
  beta_raw = "Raw beta vs generated beta",
  beta_id = "Identified beta vs generated beta",
  gamma_raw_t = "Raw Gamma vs identified truth",
  gamma_id_t = "Identified Gamma vs identified truth",
  gamma_raw = "Raw Gamma vs generated Gamma",
  gamma_id = "Identified Gamma vs generated Gamma",
  gamma_low = "Raw Gamma: least observed rows",
  gamma_high = "Raw Gamma: most observed rows"
)
curves <- rows %>%
  pivot_longer(all_of(names(panels)), names_to = "panel",
               values_to = "error") %>%
  group_by(panel, mechanism, n) %>%
  summarise(mean = mean(error), .groups = "drop") %>%
  mutate(panel = factor(panel, levels = names(panels), labels = panels))

mechanism_colors <- c("grey35", "#4fbcd8", signif_blue, "#00485b",
                      "#f9b96b", signif_orange, "#9c5200")
mechanism_shapes <- c(16, 17, 17, 17, 15, 15, 15)
legend_title <- expression("Mechanism" ~ (pi[max] / pi[min]) * ":")
# the rows on which the errors of the fitted matrix and of Gamma are measured
figure_note <- sprintf(paste(
  "Fitted matrix and Gamma: first %d rows, the rows common to every size.",
  "Beta: all its entries."
), min(curves$n))
figure <- ggplot(curves, aes(n, mean, color = mechanism, shape = mechanism)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = mechanism_colors) +
  scale_shape_manual(values = mechanism_shapes) +
  scale_x_continuous(breaks = sort(unique(curves$n))) +
  facet_wrap(~panel, scales = "free_y", ncol = 4) +
  labs(x = "Number of rows (n)", y = "RRMSE", color = legend_title,
       shape = legend_title, caption = figure_note) +
  theme_significance() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold", size = 10),
        plot.caption = element_text(hjust = 0, size = 11))

dir.create(file.path(out_dir, "paper"), showWarnings = FALSE)
ggsave(file.path(out_dir, "paper", "figure_sens_missingness.pdf"), figure,
       device = cairo_pdf, width = 14, height = 13, units = "in")

# Final figures ----
# One figure for the mechanisms driven by X and one for those driven by X
# and Z, each with MCAR as the baseline. Six panels: the fitted matrix on
# the first rows (all entries, unobserved entries), then the raw and the
# identified beta, then the raw and the identified Gamma, both against the
# generated coefficients and on all their entries. The two panels of a row
# share their axis, and the two figures have the same axes.
final_panels <- c(
  theta_all = "Fitted matrix: all entries",
  theta_unobs = "Fitted matrix: unobserved entries",
  beta_raw = "Raw beta",
  beta_id = "Identified beta",
  gamma_raw_full = "Raw Gamma",
  gamma_id_full = "Identified Gamma"
)
final_curves <- rows %>%
  pivot_longer(all_of(names(final_panels)), names_to = "panel",
               values_to = "error") %>%
  group_by(panel, mechanism, n) %>%
  summarise(mean = mean(error), .groups = "drop") %>%
  mutate(
    # the row of the panel: "theta", "beta" or "gamma"
    pair = sub("_.*", "", panel),
    panel = factor(panel, levels = names(final_panels), labels = final_panels)
  )
# the smallest and the largest mean of each row of panels, over the seven
# mechanisms: drawn as invisible points, they set the axes
final_ranges <- final_curves %>% group_by(pair) %>%
  mutate(low = min(mean), high = max(mean)) %>% ungroup() %>%
  distinct(panel, low, high) %>%
  pivot_longer(c(low, high), values_to = "mean") %>%
  mutate(n = min(final_curves$n))
final_note <- sprintf(paste(
  "Fitted matrix: first %d rows, the rows common to every size.",
  "Beta and Gamma: all their entries,\nagainst the generated coefficients."
), min(final_curves$n))

# The figure of one family: MCAR and the three mechanisms at the positions
# `cases` among the seven.
final_figure <- function(cases) {
  keep <- final_curves %>%
    filter(as.integer(mechanism) %in% c(1, cases)) %>%
    mutate(mechanism = droplevels(mechanism))
  ggplot(keep, aes(n, mean)) +
    geom_blank(data = final_ranges) +
    geom_line(aes(color = mechanism), linewidth = 1) +
    geom_point(aes(color = mechanism, shape = mechanism), size = 2.5) +
    scale_color_manual(values = mechanism_colors[c(1, cases)]) +
    scale_shape_manual(values = mechanism_shapes[c(1, cases)]) +
    scale_x_continuous(breaks = sort(unique(final_curves$n))) +
    facet_wrap(~panel, scales = "free_y", ncol = 2) +
    labs(x = "Number of rows (n)", y = "RRMSE", color = legend_title,
         shape = legend_title, caption = final_note) +
    theme_significance() +
    theme(legend.position = "top", legend.justification = "left",
          strip.text = element_text(face = "bold", size = 10),
          plot.caption = element_text(hjust = 0, size = 10,
                                      margin = margin(t = 12)))
}
ggsave(file.path(out_dir, "paper", "figure_sens_missingness_x.pdf"),
       final_figure(2:4), device = cairo_pdf, width = 8, height = 10,
       units = "in")
ggsave(file.path(out_dir, "paper", "figure_sens_missingness_xz.pdf"),
       final_figure(5:7), device = cairo_pdf, width = 8, height = 10,
       units = "in")

# Figures of the paper ----
# The figures above are for diagnosis. The paper keeps two figures, one for
# each family of mechanisms, with MCAR as the baseline in both. Three panels:
# the fitted matrix on the unobserved entries of all the rows of the fit,
# then the raw beta and the raw Gamma against the generated coefficients, on
# all their entries.
# As in identification_results.R: the labels use the symbols of the paper
# (the algorithm returns Theta hat, beta hat and Gamma hat), the theme and
# the colours are those of RSSthemes, and a figure is 8 in wide: reduced to
# the width of the text (about 6.4 in), its smallest text is about 8 pt.
# Each panel has its own axis.
theme_journal <- theme_significance() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0),
    # four keys and their title on one line of the 8 in
    legend.key.width = unit(1.8, "lines"),
    legend.title = element_text(size = 12),
    legend.text = element_text(hjust = 0, size = 11),
    strip.text = element_text(size = 13),
    panel.spacing = unit(1.1, "lines")
  )
paper_panels <- c(
  theta_unobs_full = "widehat(Theta)[bar(Omega)]",
  beta_raw = "widehat(beta)",
  gamma_raw_full = "widehat(Gamma)"
)
paper_curves <- rows %>%
  pivot_longer(all_of(names(paper_panels)), names_to = "panel",
               values_to = "error") %>%
  group_by(panel, mechanism, n) %>%
  summarise(mean = mean(error), .groups = "drop") %>%
  mutate(panel = factor(panel, levels = names(paper_panels),
                        labels = paper_panels))
# the mechanisms as the paper names them, followed by pi_max / pi_min
paper_labels <- with(mechanisms, case_when(
  kappa_x == 0 ~ "MCAR",
  kappa_z == 0 ~ sprintf("MAR(X), %.0f", pi_ratio),
  TRUE ~ sprintf("MAR(XZ), %.0f", pi_ratio)
))
# MCAR is a circle; the three strengths of a family are a triangle, a square
# and a diamond, from the weakest to the strongest (filled shapes of one size)
paper_shapes <- c(21, 24, 22, 23)

# The figure of one family: MCAR and the three mechanisms at the positions
# `cases` among the seven.
paper_figure <- function(cases) {
  keep <- paper_curves %>%
    filter(as.integer(mechanism) %in% c(1, cases)) %>%
    mutate(mechanism = droplevels(mechanism))
  ggplot(keep, aes(n, mean, color = mechanism, fill = mechanism,
                   shape = mechanism)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.4) +
    # the same four colours in the two figures: grey and three steps of blue
    scale_color_manual(values = mechanism_colors[1:4],
                       labels = paper_labels[c(1, cases)]) +
    scale_fill_manual(values = mechanism_colors[1:4],
                      labels = paper_labels[c(1, cases)]) +
    scale_shape_manual(values = paper_shapes,
                       labels = paper_labels[c(1, cases)]) +
    scale_x_continuous(breaks = sort(unique(paper_curves$n))) +
    facet_wrap(~panel, nrow = 1, scales = "free_y", labeller = label_parsed) +
    labs(x = expression(italic(n)), y = "RRMSE",
         color = legend_title, fill = legend_title, shape = legend_title) +
    theme_journal
}
ggsave(file.path(out_dir, "paper", "figure_sens_missingness_paper_x.pdf"),
       paper_figure(2:4), device = cairo_pdf, width = 8, height = 3.8,
       units = "in")
ggsave(file.path(out_dir, "paper", "figure_sens_missingness_paper_xz.pdf"),
       paper_figure(5:7), device = cairo_pdf, width = 8, height = 3.8,
       units = "in")
