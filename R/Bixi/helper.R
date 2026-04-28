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
  
    require(BKTR)
    bkdat <- BixiData$new()
    
    filename_base <- paste0(
      round(miss_pct * 100), "percent_",
      timestamp, "_",
      split_id, "_train",
      train_size
    )
    train_df <- rw_a_file(
      paste0(filename_base, "_train.rds"),
      directory = split_dir,
      type = "read"
    ) |> dplyr::rename(location = column, time = row)
    
    
    bkdat$temporal_positions_df %<>%
      filter(time %in% train_df$time)
    
    
    p_lgth <- KernelParameter$new(value = 7, is_fixed = TRUE)
    se_lgth <- KernelParameter$new(value = 6.427, is_fixed = TRUE)
    per_lgth <- KernelParameter$new(value = 1.039, is_fixed = TRUE)
    temporal_kernel <- KernelSE$new(lengthscale = se_lgth) *
      KernelPeriodic$new(lengthscale = per_lgth, period_length = p_lgth)
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
