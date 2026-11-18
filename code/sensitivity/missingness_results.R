# Sensitivity analysis, missingness mechanism: tables and figure.
#
# Reads results_missingness.rds (written by missingness_run.R), prints the
# tables and saves the figure in the folder paper/ of the run.
# Run from the root of the simulation repository:
#   Rscript code/sensitivity/missingness_results.R         (draft or full run)
#   Rscript code/sensitivity/missingness_results.R lite    (the test run)
#
# One row of the results is one fit. All errors are relative,
# ||estimate - truth||_F / ||truth||_F:
#   theta_all, theta_unobs: the fitted matrix on the first rows (the rows
#     that every size contains), on all their entries and on the unobserved
#     ones;
#   beta_raw, beta_id: the raw and the identified estimates of beta.
library(dplyr)
library(tidyr)
library(ggplot2)
library(RSSthemes)

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
errors <- c("theta_all", "theta_unobs", "beta_raw", "beta_id")
show("Errors: mean (SE)",
     rows %>% group_by(mechanism, n) %>%
       summarise(across(all_of(errors), mean_se), .groups = "drop"))

# The same errors minus their MCAR value in the same replication and at the
# same n. The SE is small because the fits are paired.
show("Difference from MCAR within a replication: mean (SE)",
     rows %>% group_by(seed, n) %>%
       mutate(across(all_of(errors), ~ .x - .x[kappa == 0])) %>%
       ungroup() %>% filter(kappa > 0) %>% group_by(mechanism, n) %>%
       summarise(across(all_of(errors), mean_se), .groups = "drop"))

# Tuning and checks ----
# lambda_beta_zero: share of the fits whose selected lambda_beta is 0.
# lambda_beta_rel: the selected value over the top of its grid.
# beta_raw_id: relative distance between the raw and the identified beta.
# beta_gap: relative distance between beta and the identified truth of the
#   first n rows, (X'X)^{-1} X' theta. It is zero at the largest size, where
#   X'A = 0. At the other sizes the identified estimate targets that truth,
#   so this distance is part of its error against beta.
show("Tuning and checks: means",
     rows %>% group_by(mechanism, n) %>%
       summarise(
         lambda_beta_zero = mean(lambda_beta == 0),
         lambda_beta_rel = mean(lambda_beta / lambda_beta_max),
         beta_raw_id = mean(beta_raw_id), lambda_m = mean(lambda_m),
         rank = mean(rank_fit), iterations = mean(n_iter),
         converged = mean(converged), seconds = mean(secs),
         beta_gap = mean(beta_gap),
         .groups = "drop"
       ))

# Figure ----
# Error against n, one line per mechanism: the fitted matrix on the first
# rows, and the identified beta. The theme and the colours are those of
# 4_plots.R: four steps of one hue for the ordered strengths.
panels <- c(theta_all = "Fitted matrix, first rows",
            beta_id = "Identified coefficients")
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
  facet_wrap(~panel, scales = "free_y") +
  labs(x = "Number of rows (n)", y = "RRMSE", color = legend_title,
       shape = legend_title) +
  theme_significance() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold", size = 11))

dir.create(file.path(out_dir, "paper"), showWarnings = FALSE)
ggsave(file.path(out_dir, "paper", "figure_sens_missingness.pdf"), figure,
       device = cairo_pdf, width = 10, height = 4.5, units = "in")
