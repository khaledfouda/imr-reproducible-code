library(dplyr)
library(tidyr)
library(IMR)
library(tidyverse)
library(magrittr)
library(kableExtra)
library(data.table)
library(bench)
library(lubridate)
library(BKTR)
#------------------------------------------------------------------------------
# configurations >>
TRAIN_SEQ <- seq(55, 75, 5)
source("./code/Bixi/config_default.R")
if (file.exists("./code/Bixi/config.R")) {
  source("./code/Bixi/config.R")
}



# the following file generates the splits. Used by 1_generate_train_test_splits.R
bixi_generate_one_split <- function(split_id = 1,
                                    file_override = FALSE,
                                    create_folder = FALSE) {
  # Load required libraries

  # define some constants -----------------------------
  miss_pct <- 0.25
  timestamp <- "Feb_last"
  seed <- 2025 + split_id
  decreasing_train <- TRUE
  train_n_steps <- 5
  train_stepsize <- 0.05
  out_dir <- "./data/Bixi/train_test_splits/"
  #----------------------------------------------
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # 1. Load raw data
  bixi_data <- BKTR::BixiData$new()
  data_df <- bixi_data$data_df
  
  # 2. Remove rows/locations with all-missing departures
  data_df <- data_df |>
    dplyr::group_by(time) |>
    dplyr::filter(!all(is.na(nb_departure))) |>
    dplyr::ungroup() |>
    dplyr::group_by(location) |>
    dplyr::filter(!all(is.na(nb_departure))) |>
    dplyr::ungroup()
  
  # 3. Select covariates & reshape for matrix input
  z_vars <- c("mean_temp_c", "total_precip_mm", "holiday", "max_temp_f", "humidity")
  
  data_df <- data_df |>
    dplyr::rename(column = location, row = time, y = nb_departure) |>
    dplyr::arrange(row, column)
  
  # Determine columns to split_id with "z_" and "x_"
  current_cols <- names(data_df)
  cols_to_z <- setdiff(current_cols, c(z_vars, "row", "column", "y"))
  
  data_df <- data_df |>
    dplyr::rename_with(\(x) paste0("z_", x), dplyr::all_of(cols_to_z)) |>
    dplyr::rename_with(\(x) paste0("x_", x), dplyr::any_of(z_vars))
  
  # 4. Discard stations with very few observations
  low_obs_stations <- c(
    "6194 - Métro Atwater (Atwater / Ste-Catherine)",
    "6019 - Métro Sherbrooke (de Rigaud / Berri)",
    "6036 - de la Commune / St-Sulpice",
    "6181 - Clark / Rachel",
    "6157 - de Brébeuf / du Mont-Royal",
    "6227 - de l'Esplanade / Laurier",
    "6136 - Métro Laurier (Rivard / Laurier)",
    "6184 - Métro Mont-Royal (Rivard / du Mont-Royal)"
  )
  
  data_df <- data_df |>
    dplyr::filter(!(column %in% low_obs_stations))
  
  # Print matrix dimensions
  num_rows <- dplyr::n_distinct(data_df$row)
  num_columns <- dplyr::n_distinct(data_df$column)
  message("Response data matrix dimension is: ", num_rows, " x ", num_columns)
  
  # 5. Initialize train/test dataframes
  train_df <- data_df |> dplyr::mutate(row_id = dplyr::row_number(), orig_y = y)
  test_df <- data_df |> dplyr::mutate(row_id = dplyr::row_number(), orig_y = y)
  
  n_total <- nrow(train_df)
  n_orig_na <- sum(is.na(train_df$orig_y))
  n_target_na <- floor(miss_pct * n_total)
  n_to_mask <- n_target_na - n_orig_na
  
  if (n_to_mask < 0) {
    stop("Already more than the target missing percentage in train_df; nothing to do.")
  }
  if (n_to_mask > (n_total - n_orig_na)) {
    stop("Not enough non-missing values to reach the target missing percentage.")
  }
  
  # 6. Sample IDs to mask
  mask_ids <- train_df |>
    dplyr::filter(!is.na(orig_y)) |>
    dplyr::slice_sample(n = n_to_mask) |>
    dplyr::pull(row_id)
  
  # Apply masking: train gets NA, test gets original values
  train_df <- train_df |>
    dplyr::mutate(y = dplyr::if_else(row_id %in% mask_ids, NA_real_, orig_y)) |>
    dplyr::select(-orig_y)
  
  test_df <- test_df |>
    dplyr::mutate(y = dplyr::if_else(row_id %in% mask_ids, orig_y, NA_real_)) |>
    dplyr::select(-orig_y)
  
  # 7. Generate sequentially decreasing training sets (if requested)
  if (decreasing_train) {
    train_seq <- round(seq(1 - miss_pct, by = -train_stepsize, length.out = train_n_steps) * 100)
    current_train <- train_df
    
    # Save the first step
    filename <- paste0("split_", split_id, "_train_", train_seq[1], ".rds")
    
    rw_a_file(filename,
              data = current_train, file_override = file_override,
              create_folder = create_folder, directory = out_dir, type = "write"
    )
    
    # Loop to mask remaining steps
    for (i in 2:length(train_seq)) {
      step_mask_count <- floor(train_stepsize * nrow(current_train))
      
      step_mask_ids <- current_train |>
        dplyr::filter(!is.na(y)) |>
        dplyr::slice_sample(n = step_mask_count) |>
        dplyr::pull(row_id)
      
      current_train <- current_train |>
        dplyr::mutate(y = dplyr::if_else(row_id %in% step_mask_ids, NA_real_, y))
      
      filename <- paste0("split_", split_id, "_train_", train_seq[i], ".rds")
      
      rw_a_file(filename,
                data = current_train, file_override = file_override,
                create_folder = create_folder, directory = out_dir, type = "write"
      )
    }
  }
  
  stopifnot("Overlap detected between train and test sets!" = sum(!is.na(train_df$y) & !is.na(test_df$y)) == 0)
  message("Percentage of observations in training set: ", round(100 * mean(!is.na(train_df$y)), 1), "%")
  message("Percentage of observations in test set: ", round(100 * mean(!is.na(test_df$y)), 1), "%")
  
  # Save main test split (main train split is unnecessary and has been removed)
  filename_base <- paste0("split_", split_id)
  
  rw_a_file(paste0(filename_base, "_test.rds"),
            data = test_df,
            file_override = file_override, create_folder = create_folder,
            directory = out_dir, type = "write"
  )
}
#-----------------------------------------------------------
# loads the saved splits and sends it the models for training. (used by fit functions)
bixi_load_split <- function(split_id = 1,
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
  #----------------------------------------------
  set.seed(seed)
  filename_base <- paste0("split_", split_id)
  
  train_df <- rw_a_file(
    paste0(filename_base, "_train_", train_size, ".rds"),
    directory = split_dir, 
    type = "read"
  ) |>
    dplyr::mutate(row = as.Date(row))
  
  test_df <- rw_a_file(
    paste0(filename_base, "_test.rds"),
           directory = split_dir,
           type = "read") |>
    dplyr::mutate(row = as.Date(row))
  
  
  
  #  Build Covariate Matrices X and Z
  X <- train_df |>
    dplyr::select(row, dplyr::starts_with("x_")) |>
    dplyr::distinct(row, .keep_all = TRUE) |>
    dplyr::select(-row) |> 
    dplyr::select(x_mean_temp_c, x_total_precip_mm)
  
  Z <- train_df |>
    dplyr::select(column, dplyr::starts_with("z_")) |>
    dplyr::distinct(column, .keep_all = TRUE) |>
    dplyr::select(-column) |> 
    dplyr::select(z_area_park)
  
  
  
  # 4. Build Response Matrix Y
  Y <- train_df |>
    dplyr::select(row, column, y) |>
    tidyr::pivot_wider(names_from = column, values_from = y) |>
    dplyr::select(-row) |>
    as.matrix()
  
  # 5. Calculate and print missing value percentages
  col_na <- colSums(!is.na(Y))
  row_na <- rowSums(!is.na(Y))
  
  message(sprintf(
    "Observed-value counts (pct):\n Columns   : min = %d (%.1f%%), max = %d (%.1f%%)\n Rows      : min = %d (%.1f%%), max = %d (%.1f%%)\n Train     : %.1f%%",
    min(col_na), 100 * min(col_na) / nrow(Y),
    max(col_na), 100 * max(col_na) / nrow(Y),
    min(row_na), 100 * min(row_na) / ncol(Y),
    max(row_na), 100 * max(row_na) / ncol(Y),
    100 * mean(!is.na(Y))
  ))
  
  # 6. Build Test Matrix
  test_set <- test_df |>
    dplyr::select(row, column, y) |>
    tidyr::pivot_wider(names_from = column, values_from = y) |>
    dplyr::select(-row) |>
    as.matrix() |>
    IMR::as_incomplete()
  
  message(sprintf("Test      : %.1f%%", 100 * sum(test_set != 0) / length(test_set)))
  
  # 7. Return final structured object
  list(
    Y = Y,
    test_mask = IMR::as_incomplete((test_set != 0) * 1),
    test = test_set,
    X = X,
    Z = Z,
    col_names = colnames(Y),
    row_names = unique(train_df$row)
  )
}
#-----------------------------------------------
# loads/generate the kernels (covariance matrices of the temporal/spatial correlation) (used by fit functions)
bixi_load_kernels <- function(split_id = 1,
                              train_size = 70,
                              return_distance = FALSE){
  
  # Validate inputs
  if (!(train_size %in% seq(55, 75, 5))) {
    stop("train_size must be one of: 55, 60, 65, 70, 75")
  }
  # define some constants -----------------------------
  miss_pct <- 0.25
  timestamp <- "Feb_last"
  seed <- 2025 + split_id
  split_dir <- "./data/Bixi/train_test_splits/"
  #----------------------------------------------
  
    bkdat <- BixiData$new()
    
    filename_base <- paste0("split_", split_id)
    
    train_df <- rw_a_file(
      paste0(filename_base, "_train_", train_size, ".rds"),
      directory = split_dir,
      type = "read"
    ) |> dplyr::rename(location = column, time = row)
    
    
    bkdat$temporal_positions_df %<>%
      filter(time %in% train_df$time)
    
    
    p_lgth <- BKTR::KernelParameter$new(value = 7, is_fixed = TRUE)
    se_lgth <- BKTR::KernelParameter$new(value = 6.427, is_fixed = TRUE)
    per_lgth <- BKTR::KernelParameter$new(value = 1.039, is_fixed = TRUE)
    temporal_kernel <- BKTR::KernelSE$new(lengthscale = se_lgth) *
      BKTR::KernelPeriodic$new(lengthscale = per_lgth, period_length = p_lgth)
    temporal_kernel$set_positions(bkdat$temporal_positions_df)
    
    bkdat$spatial_positions_df %<>%
      filter(location %in% train_df$location)
    sp_lgth <- KernelParameter$new(value = 0.018, is_fixed = TRUE)
    spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 5,lengthscale = sp_lgth)
    spatial_kernel$set_positions(bkdat$spatial_positions_df)
    
    distance = list()
    if(return_distance)
      distance = list(
        spatial = as.matrix(spatial_kernel$distance_matrix),
        temporal = as.matrix(temporal_kernel$distance_matrix)
      )
    
    temporal_kernel <- temporal_kernel$kernel_gen() %>% as.matrix()
    spatial_kernel = spatial_kernel$kernel_gen() %>% as.matrix()
    
    list(spatial=spatial_kernel, temporal=temporal_kernel, distance=distance)
  }
#-----------------------------------------------------
#  fits the BKTR model to a single data file and saves to a file. used in 2_fit_BKTR.R
bixi_fit_bktr <- function(split_id = 1,
                          train_size = 70,
                          burn_in_iter = 1000,
                          sampling_iter = 500) {
  
  # Validate inputs
  if (!(train_size %in% seq(55, 75, 5))) {
    stop("train_size must be one of: 55, 60, 65, 70, 75")
  }
  # define some constants -----------------------------
  miss_pct <- 0.25
  timestamp <- "Feb_last"
  seed <- 2025 + split_id
  split_dir <- "./data/Bixi/train_test_splits/"
  out_dir <- "./output/Bixi/model_fits/BKTR_all_fits/"
  #----------------------------------------------
  set.seed(seed)
  TSR$set_params(seed = seed, fp_type = "float32", fp_device = "cpu")
  
  # Base filename for I/O
  filename_base <- paste0("split_", split_id)
  
  # Load base objects and splits
  bixi_data <- BixiData$new()
  
  train_df <- rw_a_file(
    paste0(filename_base, "_train_", train_size, ".rds"),
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
  p_lgth <- BKTR::KernelParameter$new(value = 7, is_fixed = TRUE)
  k_local_periodic <- BKTR::KernelSE$new() * BKTR::KernelPeriodic$new(period_length = p_lgth)
  
  # Initialize BKTR model
  bktr_fit <- BKTR::BKTRRegressor$new(
    formula = y ~ 1 + x_mean_temp_c + z_area_park + x_total_precip_mm,
    data_df = train_df,
    spatial_positions_df = bixi_data$spatial_positions_df,
    temporal_positions_df = bixi_data$temporal_positions_df,
    rank = 8,
    spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 5),
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
    paste0("bktr_fit_", filename_base, "_train_", train_size, ".rds"),
    data = return_obj,
    file_override = TRUE,
    create_folder = TRUE,
    directory = out_dir,
    type = "write"
  )
  
  return(return_obj)
}

#--------------------------------------------------------
# uses the fitted BKTR model to generate the results dataframe. used in 2_fit_BKTR.R
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
  fit_dir <- "./output/Bixi/model_fits/BKTR_all_fits/"
  #----------------------------------------------
  set.seed(seed)
  
  # read train/test data
  filename_base <- paste0("split_", split_id)
  
  train_df <- rw_a_file(
    paste0(filename_base, "_train_", train_size, ".rds"),
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
  # keep this
  fit <- rw_a_file(
    paste0("bktr_fit_", filename_base, "_train_", train_size, ".rds"),
    directory = fit_dir,
    type = "read"
  )
  
  
  
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