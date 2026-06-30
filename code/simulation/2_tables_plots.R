source("./code/helper.R")
source("./code/simulation/helper.R")




results1 <- rw_a_file("results_scenario_1.rds",
                    directory = "./output/simulation/",
                    type = "read"
)
results2_1 <- rw_a_file("results_scenario_2_part_1.rds",
                                      directory = "./output/simulation/",
                                      type = "read"
)
results2_2 <- rw_a_file("results_scenario_2_part_2.rds",
                    directory = "./output/simulation/",
                    type = "read"
)

#----------------------------------------
# Generate Table 1 of Scenario 1
#----------------------------------------
compute.mean.sd <- function(x){
  s <- sd(x)
  if(is.na(s)) "" else paste0(round(mean(x),3)," (", round(s,3) ,")")
}

results1 %>%
  dplyr::select(train_rrmse, test_rrmse, beta_rrmse, 
                M_rrmse,
                rank, model, dim) %>%
  arrange(dim, model) %>%
  dplyr::group_by(dim, model) %>%
  summarize_all(compute.mean.sd)  %>%
  rename(beta = beta_rrmse, 
         M = M_rrmse,
         train = train_rrmse,
         test = test_rrmse) %>%
  pivot_longer(c("beta", "M",  "train", "test", "rank"),
               names_to = c("metric")) %>%
  mutate(value = if_else(model == "SI" & metric == "M","",value)) %>%
  mutate(model = if_else(model == "SI", "SoftImpute", model)) ->
  sim1.tab



metric_labels <- c(
  beta = "RRMSE($\\beta$)",
  M    = "RRMSE($M$)",
  #theta    = "RMSE($\\Theta$)",
  train = "RRMSE(train)",
  test = "RRMSE(test)",
  rank       = "Rank"
)

method_order <- c(
  "IMR","SoftImpute",  "MCCI"
)

wide <- sim1.tab %>%
  mutate(metric = recode(metric, !!!metric_labels)) %>%
  mutate(model = factor(model, levels = method_order)) %>%
  arrange(dim, model) %>%
  dplyr::select(dim, model, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  as.data.frame()

# How many rows per panel
panel_counts <- wide %>% count(dim, name = "rows")
panel_index  <- stats::setNames(panel_counts$rows,
                                paste0("n = m = ", panel_counts$dim))
rows_to_bold <- c(1, 4, 7, 10)
bold_vector <- seq_len(nrow(wide)) %in% rows_to_bold
# Build table
kbl(
  wide %>% dplyr::select(-dim),
  format    = "latex",
  booktabs  = TRUE,
  escape    = FALSE,
  label =  "sim1",
  position  = "!tb",
  col.names = c(
    "Model", 
    "RRMSE($\\beta$)", 
    "RRMSE($M$)", 
    "RRMSE($\\btheta_{\\text{train}}$)", 
    "RRMSE($\\btheta_\\text{test}$)", 
    "Rank($\\btheta$)"
  ),
  align     = c("l", rep("r", ncol(wide) - 2)),
  caption   = paste(
    "Empirical relative root mean squared errors (RRMSEs),",
    "estimated ranks, and standard deviations (in parentheses)",
    "under the model $\\btheta=\\bX\\bbeta+\\bM$.",
    "We consider dimensions $(n,m)=(400,400),(600,600),(800,800),(1000,1000)$",
    "and a fixed missingness rate of $80\\%$. The true rank of $\\btheta$ is 14."
  )) |>
  pack_rows(index = panel_index, bold = FALSE) |>
  kable_styling(latex_options = c("hold_position")) |>
  column_spec(5, color = "black",bold = bold_vector) |>
  column_spec(1, color = "black",bold = TRUE ) |>
  row_spec(1:nrow(wide), color = "black") -> sim1.tbl;sim1.tbl
#-------------------------------------------------------------------------------
# Generate Figure 1 of Scenario 2
#-------------------------------------------------------------------------------
plot_data <- results2_1 |>
  mutate(
    sparsity = round(missing_pct, 2),
    model = if_else(model == "SI", "Soft-Impute", model),
    M_rrmse = if_else(model == "Soft-Impute", NA_real_, M_rrmse)
  ) |>
  select(
    theta = theta_rrmse, test = test_rrmse, beta = beta_rrmse, 
    M = M_rrmse, gamma = gamma_rrmse, rank, sparsity, model
  ) |>
  mutate(rank = rank - 15) |> 
  summarise(
    across(
      c(theta, test, beta, M, gamma, rank), 
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))
    ),
    .by = c(model, sparsity)
  ) |>
  pivot_longer(
    cols = -c(model, sparsity),
    names_to = c("metric", ".value"),
    names_sep = "_" 
  ) |>
  mutate(
    ymin = pmax(0, mean - sd), 
    ymax = mean + sd
  )
metric_labels <- c(
  beta  = "RRMSE(beta)",
  gamma = "RRMSE(Gamma)",
  M     = "RRMSE(M)",
  theta = "RRMSE(Theta)",
  test  = "RRMSE(Theta[plain(test)])",
  rank  = "textstyle('Estimated Rank') - textstyle('True Rank')"
)

plot_data <- plot_data |>
  mutate(
    metric_lab = factor(metric, levels = names(metric_labels), labels = unname(metric_labels))
  )

bounds_data <- data.frame(
  metric_lab = factor(
    rep(metric_labels[1:5],each=2),
    levels = unname(metric_labels)
  ),
  mean = c(0.4, 1.0, 0.4, 1.0, 0.4, 1.0, 0.2, 0.7, 0.2, 0.7),
  sparsity = 0.7 
)
rank_line_data <- data.frame(
  metric_lab = factor(metric_labels[6], levels = unname(metric_labels)),
  yintercept = 15
)
pcolors <- c(signif_blue, signif_orange)

sim2.g <- ggplot(plot_data, aes(x = sparsity, y = mean, color = model, fill = model)) +
    geom_blank(data = bounds_data, aes(y = mean), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  
  scale_color_manual(values = pcolors) +
  scale_fill_manual(values  = pcolors) +
  scale_y_continuous(n.breaks = 5) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks =seq(0.7, 0.95, by = 0.05)) +
  facet_wrap(~ metric_lab, labeller = label_parsed, ncol = 3, scales = "free_y") +
  labs(x = "Sparsity Level", y = "Estimation Error", color = "Model", fill = "Model") +
  theme_significance() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    strip.text = element_text(face = "bold", size = 11)
  ); sim2.g



ggsave(
  filename = "./output/simulation/figure_1_scenario_2.pdf", 
  plot = sim2.g, 
  device = cairo_pdf,  
  width = 12,           
  height = 6,          
  units = "in"
)
#-------------------------------------------------------------------------------
# Generate Figure 2 of Scenario 2
#-------------------------------------------------------------------------------



plot_titles <- c("Error Ratio (IMR / Soft-Impute)", "Time Ratio (IMR / Soft-Impute)")

plot_data <- results2_2 |>
  mutate(sparsity = round(missing_pct, 2)) |> 
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    mean_error = mean(test_rrmse, na.rm = TRUE),
    .by = c(model, sparsity)
  ) |>
  pivot_wider(names_from = model, values_from = c(mean_time, mean_error)) |>
  transmute(
    sparsity,
    `Time Ratio (IMR / Soft-Impute)` = mean_time_IMR / mean_time_SI,
    `Error Ratio (IMR / Soft-Impute)` = mean_error_IMR / mean_error_SI
  ) |>
  pivot_longer(
    cols = -sparsity,
    names_to = "measure_label",
    values_to = "value"
  ) |>
  mutate(measure_label = factor(measure_label, levels = plot_titles)) # to have the order i want


limit_time <- c(0, 0.4)
limit_error <- c(0, 1.0)

facet_bounds <- data.frame(
  measure_label = rep(plot_titles, each=2),
  sparsity = 0.7,            
  value = c(limit_error, limit_time) 
)

diff_plot <- ggplot(plot_data, aes(x = sparsity, y = value)) +
  geom_blank(data = facet_bounds) +
  geom_hline(aes(yintercept = 1), data.frame(measure_label=plot_titles[1]), linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1, color=signif_blue) +
  geom_point(size = 2) +
  facet_wrap(~ measure_label, scales = "free_y", ncol = 1) +
  scale_x_continuous(
    name = "Sparsity Level", 
    labels = percent_format(accuracy = 1),
    breaks = seq(0.7, 0.95, by = 0.05)
  ) +
  scale_y_continuous(
    name = "Ratio",
    breaks = ~ if (max(.x, na.rm = TRUE) <= 0.4) {
      seq(limit_time[1], limit_time[2], by = 0.2) 
    } else {
      seq(limit_error[1], limit_error[2], by = 0.2) # Breaks for the Error Ratio 
    }
  ) +
  theme_significance() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95")
  ); diff_plot


ggsave(
  filename = "./output/simulation/figure_2_scenario_2.pdf", 
  plot = diff_plot, 
  device = cairo_pdf,  
  width = 12,           
  height = 6,          
  units = "in"
)

#----------------------------------
# DONE

