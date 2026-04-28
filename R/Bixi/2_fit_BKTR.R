source("./R/helper.R")
source("./R/Bixi/helper.R")

bixi_fit_bktr <- function(split_id = 1,
                          train_size = 70,
                          burn_in_iter = 1000,
                          sampling_iter = 500) {
  require(BKTR)
  require(data.table)
  require(dplyr)
  require(bench)
  require(lubridate)

  # Validate inputs
  if (!(train_size %in% seq(55, 75, 5))) {
    stop("train_size must be one of: 55, 60, 65, 70, 75")
  }
  # define some constants -----------------------------
  miss_pct <- 0.25
  timestamp <- "Feb_last"
  seed <- 2025 + split_id
  split_dir <- "./data/Bixi/train_test_splits/"
  out_dir <- "./data/Bixi/model_fits/BKTR_all_fits/"
  #----------------------------------------------
  set.seed(seed)
  TSR$set_params(seed = seed, fp_type = "float32", fp_device = "cpu")

  # Base filename for I/O
  filename_base <- paste0(
    round(miss_pct * 100), "percent_",
    timestamp, "_",
    split_id, "_train",
    train_size
  )

  # Load base objects and splits
  bixi_data <- BixiData$new()

  train_df <- rw_a_file(
    paste0(filename_base, "_train.rds"),
    directory = split_dir,
    type = "read"
  )

  # Format training data to data.table
  train_df <- train_df |> dplyr::rename(location = column, time = row)
  data.table::setDT(train_df)

  # Filter BKTR positions to match train.df
  bixi_data$temporal_positions_df <- bixi_data$temporal_positions_df |>
    dplyr::filter(time %in% train_df$time)

  bixi_data$spatial_positions_df <- bixi_data$spatial_positions_df |>
    dplyr::filter(location %in% train_df$location)

  train_df[, time := as.character(time)]
  bixi_data$temporal_positions_df[, time := as.character(time)]

  setkey(train_df, location, time)
  setkey(bixi_data$temporal_positions_df, time)
  setkey(bixi_data$spatial_positions_df, location)

  # Set up BKTR kernels
  p_lgth <- KernelParameter$new(value = 7, is_fixed = TRUE)
  k_local_periodic <- KernelSE$new() * KernelPeriodic$new(period_length = p_lgth)

  # Initialize BKTR model
  bktr_fit <- BKTRRegressor$new(
    formula = y ~ 1 + x_mean_temp_c + z_area_park + x_total_precip_mm,
    data_df = train_df,
    spatial_positions_df = bixi_data$spatial_positions_df,
    temporal_positions_df = bixi_data$temporal_positions_df,
    rank = 8,
    spatial_kernel = KernelMatern$new(smoothness_factor = 5),
    temporal_kernel = k_local_periodic,
    burn_in_iter = burn_in_iter,
    sampling_iter = sampling_iter
  )

  # Run MCMC sampling and track time
  start_time <- Sys.time()
  mcmc_bench <- bench::bench_time(bktr_fit$mcmc_sampling())

  total_time_secs <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))

  # Compile return object
  return_obj <- list(
    fit = bktr_fit,
    time = total_time_secs,
    time1 = lubridate::time_length(mcmc_bench[1], "seconds"),
    time2 = lubridate::time_length(mcmc_bench[2], "seconds")
  )

  # Save the fit
  rw_a_file(
    paste0("bktr_fit_", filename_base, ".rds"),
    data = return_obj,
    file_override = TRUE,
    create_folder = TRUE,
    directory = out_dir,
    type = "write"
  )

  return(return_obj)
}
#-------------------------------------------------------
# fit all 50 * 5 data files.
train_seq <- seq(55, 75, 5)
num_splits <- 50

for (split_id in 1:num_splits) {
  for (train_size in train_seq) {
    bktr_fit <- bixi_fit_bktr(
      split_id, train_size,
      burn_in_iter = 1000,
      sampling_iter = 500
    )
  }
}
#--------------------------------
# since we save each one of the 250 fits to a file. we combine them to a single file
# though we only combine the diagnostic output
bixi_fit_bktr_post <- function(split_id = 1,
                               train_size = 70) {
  # Validate inputs
  if (!(train_size %in% seq(55, 75, 5))) {
    stop("train_size must be one of: 55, 60, 65, 70, 75")
  }
  # define some constants -----------------------------
  miss_pct <- 0.25
  timestamp <- "Feb_last"
  seed <- 2025 + split_id
  split_dir <- "./data/Bixi/train_test_splits/"
  fit_dir <- "./data/Bixi/model_fits/BKTR_all_fits/"
  #----------------------------------------------
  set.seed(seed)

  # read train/test data
  filename_base <- paste0(
    round(miss_pct * 100), "percent_",
    timestamp, "_",
    split_id
  )
  train_df <- rw_a_file(
    paste0(filename_base, "_train", train_size, "_train.rds"),
    directory = split_dir,
    type = "read"
  ) |>
    dplyr::rename(location = column, time = row) |>
    dplyr::mutate(time = as.Date(time))

  test_df <- rw_a_file(
    paste0(filename_base, "_test.rds"),
    directory = split_dir,
    type = "read"
  ) |>
    dplyr::rename(location = column, time = row) |>
    dplyr::mutate(time = as.Date(time))

  # read BKTR fit file
  filename_base <- paste0(
    round(miss_pct * 100), "percent_",
    timestamp, "_",
    split_id, "_train",
    train_size
  )

  # keep this
  rw_a_file(
    paste0("bktr_fit_", filename_base, "_train", train_size, ".rds"),
    directory = fit_dir,
    type = "read"
  )

  # remove this later
  # filename_base <- paste0(
  #   round(miss_pct * 100), "percent_",
  #   timestamp, "_",
  #   split_id
  # )
  # fit <- rw_a_file(
  #   paste0(
  #     "bktrfit_", filename_base, "_train", split_id,
  #     "_train", train_size, "_.rds"
  #   ),
  #   directory = fit_dir,
  #   type = "read"
  # )
  #-- end of remove later

  # obtain test estimates
  fit$fit$imputed_y_estimates |>
    as.data.frame() |>
    merge(test_df, by = c("location", "time")) |>
    dplyr::select(location, time, y_est, y) ->
  test.estimates

  # obtain train estimates
  fit$fit$imputed_y_estimates |>
    as.data.frame() |>
    merge(filter(train_df, !is.na(y)), by = c("location", "time")) |>
    dplyr::select(location, time, y_est, y) ->
  train.estimates


  evaluate_estimates(train.estimates$y, train.estimates$y_est,
    test.estimates$y, test.estimates$y_est,
    time = fit$time,
    model = "BKTR"
  ) %>%
    mutate(rank_m = fit$fit$rank_decomp)
}

#--------------------------------------------------------------
train_seq <- seq(55, 75, 5)
num_splits <- 50
all_results <- data.frame()

for (split_id in 1:num_splits) {
  for (train_size in train_seq) {
    all_results <- rbind(
      all_results,
      bixi_fit_bktr_post(split_id, train_size) %>%
        mutate(
          split_id = split_id,
          train_size = train_size
        )
    )
  }
  print(paste0("split #", split_id))
}

# save to disk >> 
rw_a_file(
  "results_bktr.rds",
  data = all_results,
  file_override = TRUE,
  create_folder = TRUE,
  directory = "./data/Bixi/",
  type = "write"
)

# checking the results > 
all_results %>%
  select(-model) %>%
  group_by(train_size) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  as.data.frame()

#--- DONE

