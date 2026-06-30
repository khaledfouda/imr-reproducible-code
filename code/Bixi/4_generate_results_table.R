source("./code/helper.R")
source("./code/Bixi/helper.R")

rw_a_file(
  "results_bktr.rds",
  directory = "./output/Bixi/",
  type = "read"
) |>
  rbind(
    rw_a_file(
      "results_imr.rds",
      directory = "./output/Bixi/",
      type = "read"
    ) 
  ) |>
  group_by(model, train_size) |>
  summarise(across(everything(), list(mean = mean, sd = sd)), .groups = "drop") |>
  arrange(train_size, test_rrmse_mean) |>
  dplyr::select(model, train_size, contains("rrmse"), contains("time")) |>
  as.data.frame() |>
  dplyr::select(model, train_size, test_rrmse_mean, test_rrmse_sd, time_mean, time_sd) |>
  group_by(train_size) |>
  mutate(
    is_min_rrmse = test_rrmse_mean == min(test_rrmse_mean, na.rm = TRUE),
    is_min_time = time_mean != max(time_mean, na.rm = TRUE)
  ) -> df


df |> 
  ungroup() |>
  mutate(
    train_size = sprintf("\\textbf{%s\\%%}", train_size),
    rrmse_str = sprintf("%.4f (%.4f)", test_rrmse_mean, test_rrmse_sd),
    rrmse_str = if_else(is_min_rrmse, sprintf("\\textbf{%s}", rrmse_str), rrmse_str),
    time_str = sprintf("%.2f (%.2f)", time_mean, time_sd),
    time_str = if_else(is_min_time, sprintf("\\textbf{%s}", time_str), time_str)
  ) |>
  dplyr::select(train_size, model, rrmse_str, time_str) |>
  mutate(model = factor(model, levels = c("BKTR", "IMR-S", "IMR-N"))) |>
  pivot_wider(
    names_from = model, 
    values_from = c(rrmse_str, time_str)
  ) ->
  table_data


latex_code <-
  table_data |>
  kbl(
    format = "latex",
    escape = FALSE,       
    booktabs = TRUE,
    align = "lcccccc",     
    col.names = c(
      "{\\bfseries Train Size}",
      "{\\bfseries \\texttt{BKTR}}", 
      "{\\bfseries \\texttt{IMR-S}}", 
      "{\\bfseries \\texttt{IMR-N}}",
      "{\\bfseries \\texttt{BKTR}}", 
      "{\\bfseries \\texttt{IMR-S}}", 
      "{\\bfseries \\texttt{IMR-N}}"
    ),
    caption = paste("Predictive performance and computational efficiency on the BIXI dataset.",
      "Results are averaged over 50 independent train/test splits,",
      "with the test fraction fixed at $14\\%$, and are reported as mean (standard deviation)."),
    label = "bixi:res"
  ) |>
  kable_styling(
    position = "center",
    font_size = 9 # Natively approximates \small
  ) |>
  add_header_above(
    c(" " = 1, 
      "{\\\\bfseries Test RRMSE}" = 3, 
      "{\\\\bfseries Computation Time (Seconds)}" = 3),
    escape = FALSE
  ) |>
  row_spec(0, extra_latex_after = "\\addlinespace[2.5pt]"); print(latex_code)

#---- END

# EXTRA:
# speed difference tables for diagnostics / reporting in the manuscript.
df |>
  dplyr::group_by(train_size) |>
  dplyr::mutate(bktr_time = time_mean[model == "BKTR"]) |>
  dplyr::ungroup() |>
  dplyr::filter(model %in% c("IMR-S", "IMR-N")) |>
  dplyr::mutate(times_faster =  bktr_time / time_mean) |>
  dplyr::select(train_size, model, times_faster) %T>%
  print() |> 
  #group_by(model) |> 
  dplyr::summarise(avg_times_faster = mean(times_faster, na.rm = TRUE)) |>
  dplyr::pull(avg_times_faster)


# performance difference
df |>
  dplyr::group_by(train_size) |>
  dplyr::summarise(
    bktr_rrmse = test_rrmse_mean[model == "BKTR"],
    imrs_rrmse = test_rrmse_mean[model == "IMR-S"],
    .groups = "drop"
  ) |>
  dplyr::mutate(
    improvement_pct = (bktr_rrmse - imrs_rrmse) / bktr_rrmse * 100
  ) |>
  dplyr::select(train_size, improvement_pct)
#---------------
