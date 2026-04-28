library(dplyr)
library(tidyr)
library(kableExtra)
source("./R/helper.R")

df <- rw_a_file(
  "results_bktr.rds",
  directory = "./data/Bixi/",
  type = "read"
) |>
  rbind(
    rw_a_file(
      "results_imr.rds",
      directory = "./data/Bixi/",
      type = "read"
    )
  ) |>
  group_by(model, train_size) |>
  summarise(across(everything(), list(mean = mean, sd = sd)), .groups = "drop") |>
  arrange(train_size, test_rrmse_mean) |>
  select(model, train_size, contains("rrmse"), contains("time")) |>
  as.data.frame() |>
  select(model, train_size, test_rrmse_mean, test_rrmse_sd, time_mean, time_sd) |>
  group_by(train_size) |>
  mutate(
    is_min_rrmse = test_rrmse_mean == min(test_rrmse_mean, na.rm = TRUE),
    is_min_time = time_mean == min(time_mean, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(
    train_size = sprintf("\\textbf{%s\\%%}", train_size),
    rrmse_str = sprintf("%.4f (%.4f)", test_rrmse_mean, test_rrmse_sd),
    rrmse_str = if_else(is_min_rrmse, sprintf("\\textbf{%s}", rrmse_str), rrmse_str),
    time_str = sprintf("%.2f (%.2f)", time_mean, time_sd),
    time_str = if_else(is_min_time, sprintf("\\textbf{%s}", time_str), time_str)
  ) |>
  select(train_size, model, rrmse_str, time_str) |>
  mutate(model = factor(model, levels = c("BKTR", "IMR-S", "IMR-N"))) |>
  pivot_wider(
    names_from = model, 
    values_from = c(rrmse_str, time_str)
  ) |>
  select(
    `Train Size` = train_size,
    BKTR = rrmse_str_BKTR,
    `IMR-S` = `rrmse_str_IMR-S`,
    `IMR-N` = `rrmse_str_IMR-N`,
    `BKTR ` = time_str_BKTR,    
    `IMR-S ` = `time_str_IMR-S`,
    `IMR-N ` = `time_str_IMR-N`
  ) ->
  table_data

latex_code <- table_data |>
  kbl(
    format = "latex",
    escape = FALSE,       
    booktabs = TRUE,      
    align = "lcccccc"     
  ) |>
  kable_styling(position = "center") |>
  add_header_above(
    c(" " = 1, "Test RRMSE" = 3, "Computation Time (Seconds)" = 3),
    bold = TRUE
  ) |>
  row_spec(0, bold = TRUE) |> 
  row_spec(1:(nrow(table_data) - 1), hline_after = TRUE) 

cat(latex_code)