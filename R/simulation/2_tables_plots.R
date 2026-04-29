source("./R/helper.R")
source("./R/other_models/SoftImpute_cv.R")
source("./R/other_models/MCCI.R")
source("./R/simulation/helper.R")

require(tidyverse)
require(IMR)
require(kableExtra)
require(magrittr)
require(scales)

results1 <- rw_a_file("results_scenario_1.rds",
                    directory = "./data/Simulation/",
                    type = "read"
)
results2_1 <- rw_a_file("results_scenario_2_part_1.rds",
                                      directory = "./data/Simulation/",
                                      type = "read"
)
results2_2 <- rw_a_file("results_scenario_2_part_2.rds",
                    directory = "./data/Simulation/",
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
  select(dim, model, metric, value) %>%
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
  format    = "html",
  booktabs  = TRUE,
  escape    = FALSE,
  label =  "tab:sim1",
  position  = "htbp",
  # col.names = c("Method", setdiff(names(wide), c("n", "method"))),
  align     = c("l", rep("r", ncol(wide) - 2)),
  caption   = paste(
    "Empirical relative root mean square errors (RRMSEs),",
    "estimated ranks, and their standard deviations (in parentheses)",
    "under model $\\Theta=X\\beta+M$ and with dimensions $(n,m)=(400,400),(600,600),(800,800),(1000,1000)$,",
    "$80\\%$ missingness rate, and the true rank is 14.")
) |>
  pack_rows(index = panel_index, bold = FALSE) |>
  kable_styling(latex_options = c("hold_position", "striped"), font_size = 12) |>
  column_spec(5, color = "black",bold = bold_vector) |>
  column_spec(1, color = "black",bold = TRUE ) |>
  row_spec(1:nrow(wide), color = "black") -> sim1.tbl;sim1.tbl
#-------------------------------------------------------------------------------
# Generate Figure 1 of Scenario 2
#-------------------------------------------------------------------------------


results2_1 %>%
  mutate(sparsity = round(missing_pct, 2)) %>%
  dplyr::select(theta_rrmse, test_rrmse, beta_rrmse, 
                M_rrmse, gamma_rrmse, sparsity,
                rank, model) %>%
  rename(beta = beta_rrmse, 
         M = M_rrmse,
         theta = theta_rrmse,
         gamma = gamma_rrmse,
         test = test_rrmse) %>%
  dplyr::group_by(model, sparsity) %>%
  dplyr::mutate(M = if_else(model == "SI", NA, M)) %>%
  dplyr::mutate(model = if_else(model == "SI", "SoftImpute", model)) %>%
  dplyr::summarize_all(c(error_mean=mean,error_sd= sd)) %>%
  dplyr::ungroup() %>%
  pivot_longer(-c(model, sparsity),
               names_to = c("metric", "stat"),
               names_pattern = "^(.*)_error_(mean|sd)$",
               values_to = "val") %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  mutate(sparsity = sparsity) %>%
  arrange(model, sparsity) %>%
  mutate(ymin = pmax(0, mean-sd), ymax = mean+sd)  ->
  sim2.long


metric_labels <- c(
  beta = "RRMSE(beta)",
  gamma = "RRMSE(Gamma)",
  M    = "RRMSE(M)",
  theta    = "RRMSE(Theta)",
  test = "RRMSE(test)",
  rank       = "Estimated~Rank"
)
sim2.long %<>%
  mutate(metric_lab = factor(metric,
                             levels=names(metric_labels),
                             labels = unname(metric_labels)))

rank_line_data <- data.frame(
  metric_lab = factor("Estimated~Rank", levels = unname(metric_labels)),
  yintercept = 15
)
okabe_ito <- c("#56B4E9","#E69F00")
metrics <- unique(sim2.long$metric)

ggplot(sim2.long, aes(x = sparsity, y = mean, color = model, fill = model, group = model)) +
  geom_hline(data = rank_line_data, aes(yintercept = yintercept),
             linetype = "dashed", color = "black", alpha = 0.6) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  geom_point(size = 1.6) +
  scale_color_manual(values = okabe_ito) +
  scale_fill_manual(values  = okabe_ito) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(.7, 0.95, length=6)) +
  facet_wrap(~ metric_lab, labeller = label_parsed, ncol=3, scales="free_y") +
  labs(x = "Sparsity Level", y = NULL, color = "Model", fill = "Model") +
  theme_minimal(base_size = 10) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    strip.text = element_text(face = "bold")
  ) -> sim2.g; sim2.g

ggsave("./data/Simulation/figure_1_scenario_2.png", 
       sim2.g, width = 320/25.4, height = 150/25.4, dpi = 600)

#-------------------------------------------------------------------------------
# Generate Figure 2 of Scenario 2
#-------------------------------------------------------------------------------

results2_2 %>%
  mutate(sparsity = round(missing_pct, 2)) %>%
  dplyr::select(test_rrmse, sparsity, time,
                model) %>%
  rename(test = test_rrmse) %>%
  group_by(model, sparsity) %>%
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    mean_error = mean(test, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(sparsity, mean_time) %>%
  pivot_wider(
    names_from = model,
    values_from = c(mean_time, mean_error)
  ) %>%
  mutate(
    time_ratio = mean_time_IMR / mean_time_SI,
    error_improve = (mean_error_SI - mean_error_IMR) / mean_error_SI
  ) %>%
  select(sparsity, time_ratio, error_improve) %>%
  pivot_longer(
    cols = c(time_ratio, error_improve),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    measure_label = case_when(
      measure == "time_ratio" ~ "Relative Computational Cost",
      measure == "error_improve" ~ "RRMSE Reduction (%) Relative to SoftImpute"
    )
  ) %>%
  mutate(measure_label = factor(measure_label, levels = c(
    "RRMSE Reduction (%) Relative to SoftImpute",
    "Relative Computational Cost"
  )))-> plot_data

bounds <- data.frame(
  measure_label = factor("RRMSE Reduction (%) Relative to SoftImpute",
                         levels = levels(plot_data$measure_label)),
  sparsity = 0.95,
  value = c(0, 0.7)
)

ggplot(plot_data, aes(x = (sparsity), y = value)) +
  geom_blank(data = bounds) +
  
  geom_hline(data = filter(plot_data, measure == "time_ratio"),
             aes(yintercept = 1), linetype = "dashed", color = "grey50") +
  geom_hline(data = filter(plot_data, measure == "error_improve"),
             aes(yintercept = 0), linetype = "dashed", color = "grey50") +
  
  geom_line(size = 1, color = "#E69F00") +
  geom_point(size = 2, color = "#E69F00") +
  
  facet_wrap(~ measure_label, scales = "free_y", ncol = 1) +
  scale_x_continuous("Sparsity Level", labels = percent_format(accuracy = 1),
                     breaks = seq(0.95, 0.7, length.out=6)) +
  scale_y_continuous(
    name = NULL,
    labels = function(x) {
      ifelse(x <= 1 & x > -1, percent(x, accuracy = 1), number(x, accuracy = 0.1, suffix = "x"))
    },
    breaks = function(limits) {
      if (limits[2] <= 0.8) {
        return(c(0, .1, .3, .5, .7))#c(seq(0, 0.7, .2),.7))
      } else {
        return(scales::extended_breaks()(limits))
      }
    }
  ) +
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95")
  ) -> diff_plot; diff_plot

ggsave("./data/Simulation/figure_2_scenario_2.png",  
       diff_plot, width = 320/25.4, height = 150/25.4, dpi = 600)

#----------------------------------
# DONE

