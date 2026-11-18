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
#   beta_raw_t, beta_id_t: the raw and the identified estimates of beta,
#     against the identified truth of the rows of the fit, T(theta). It is
#     the generated beta at the largest size only;
#   beta_raw, beta_id: the same estimates against the generated beta;
#   gamma_raw_t, gamma_id_t, gamma_raw, gamma_id: the same four errors for
#     Gamma, on the first rows;
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
     rows %>% count(kappa, status, name = "fits"))

# The comparisons are paired: only the replications whose fits all ran.
rows <- rows %>% group_by(seed) %>% filter(all(status == "ok")) %>% ungroup()

# Mechanisms ----
# pi_ratio: the largest probability of observation over the smallest, on
# all the rows. kappa = 0 is MCAR.
mechanisms <- rows %>% filter(n == max(n)) %>% group_by(kappa) %>%
  summarise(
    pi_ratio = mean(pi_ratio), obs_rate = mean(obs_rate),
    fewest_in_a_row = min(min_row_obs), replications = n_distinct(seed)
  ) %>%
  mutate(mechanism = ifelse(kappa == 0, "MCAR",
                            sprintf("MAR (%.0f)", pi_ratio)))
show("Mechanisms, on the largest size", mechanisms)

rows <- rows %>%
  mutate(mechanism = factor(kappa, levels = mechanisms$kappa,
                            labels = mechanisms$mechanism))

# Errors ----
# Two tables for each group of errors: the errors, then the errors minus
# their MCAR value in the same replication and at the same n. The SE of a
# difference is small because the fits are paired.
error_groups <- list(
  "Fitted matrix" = c("theta_all", "theta_unobs", "theta_low", "theta_high"),
  "Beta" = c("beta_raw_t", "beta_id_t", "beta_raw", "beta_id"),
  "Gamma" = c("gamma_raw_t", "gamma_id_t", "gamma_raw", "gamma_id",
              "gamma_low", "gamma_high")
)
for (group in names(error_groups)) {
  errors <- error_groups[[group]]
  show(paste0(group, ", errors: mean (SE)"),
       rows %>% group_by(mechanism, n) %>%
         summarise(across(all_of(errors), mean_se), .groups = "drop"))
  show(paste0(group, ", difference from MCAR within a replication: mean (SE)"),
       rows %>% group_by(seed, n) %>%
         mutate(across(all_of(errors), ~ .x - .x[kappa == 0])) %>%
         ungroup() %>% filter(kappa > 0) %>% group_by(mechanism, n) %>%
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
# The panels of the paper are chosen from these later. The theme and the
# colours are those of 4_plots.R: four steps of one hue for the ordered
# strengths.
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

strength_colors <- c("#4fbcd8", signif_blue, "#00718e", "#00485b")
legend_title <- expression("Mechanism" ~ (pi[max] / pi[min]) * ":")
figure <- ggplot(curves, aes(n, mean, color = mechanism, shape = mechanism)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = strength_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_x_continuous(breaks = sort(unique(curves$n))) +
  facet_wrap(~panel, scales = "free_y", ncol = 4) +
  labs(x = "Number of rows (n)", y = "RRMSE", color = legend_title,
       shape = legend_title) +
  theme_significance() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold", size = 10))

dir.create(file.path(out_dir, "paper"), showWarnings = FALSE)
ggsave(file.path(out_dir, "paper", "figure_sens_missingness.pdf"), figure,
       device = cairo_pdf, width = 14, height = 13, units = "in")
