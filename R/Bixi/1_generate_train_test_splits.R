source("./R/Bixi/helper.R")

bixi_generate_one_split <- function(split_id = 1,
                                    file_override = FALSE,
                                    create_folder = FALSE) {
  # Load required libraries
  require(BKTR)
  require(dplyr)
  require(tidyr)
  
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
    new_split_id <- paste0(split_id, "_train", train_seq[1])
    filename <- paste0(round(miss_pct * 100), "percent_", timestamp, "_", new_split_id, "_train.rds")

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

      new_split_id <- paste0(split_id, "_train", train_seq[i])
      filename <- paste0(round(miss_pct * 100), "percent_", timestamp, "_", new_split_id, "_train.rds")

      rw_a_file(filename,
        data = current_train, file_override = file_override,
        create_folder = create_folder, directory = out_dir, type = "write"
      )
    }
  }

  stopifnot("Overlap detected between train and test sets!" = sum(!is.na(train_df$y) & !is.na(test_df$y)) == 0)
  message("Percentage of observations in training set: ", round(100 * mean(!is.na(train_df$y)), 1), "%")
  message("Percentage of observations in test set: ", round(100 * mean(!is.na(test_df$y)), 1), "%")

  # Save main train and test splits
  filename_base <- paste0(round(miss_pct * 100), "percent_", timestamp, "_", split_id, "_")

  rw_a_file(paste0(filename_base, "train.rds"),
    data = train_df,
    file_override = file_override, create_folder = create_folder,
    directory = out_dir, type = "write"
  )

  rw_a_file(paste0(filename_base, "test.rds"),
    data = test_df,
    file_override = file_override, create_folder = create_folder,
    directory = out_dir, type = "write"
  )
}
#----------------------------------------------------------------
# change file_override to TRUE to create new files

num_splits <- 50

for(split_id in 1:num_splits){
  bixi_generate_one_split(split_id,
                          file_override = FALSE,
                          create_folder = FALSE)
}
# DONE
#-------------------------------------------------------------------
