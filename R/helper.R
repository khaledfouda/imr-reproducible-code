rw_a_file <- function(filename,
                              data = NULL,
                              file_override = FALSE,
                              create_folder = FALSE,
                              directory = "./data/Bixi/",
                              type = "read") {
  stopifnot(type %in% c("read", "write"))
  # Handle directory creation
  if (!dir.exists(directory)) {
    if (create_folder) {
      dir.create(directory, recursive = TRUE)
    } else {
      stop("Folder does not exist: ", directory)
    }
  }
  full_path <- file.path(directory, filename)
  
  # Handle Reading
  if (type == "read") {
    if (file.exists(full_path)) {
      return(readRDS(full_path))
    } else {
      stop("Attempting to read a non-existing file: ", full_path)
    }
  } else if (type == "write") {
    # Handle Writing
    stopifnot(!is.null(data))
    if (file.exists(full_path)) {
      if (!file_override) {
        stop("File already exists: ", full_path)
      }
      message("Warning: Overwriting existing file: ", full_path)
    }
    
    saveRDS(data, file = full_path)
    invisible(full_path)
  }
}
#-----------------
comp_rank <- function(x){
  if(100 < min(dim(x)/2)){
    return(sum(irlba::irlba(x, nv=100)$d > 1e-5 ))
  }else
    return(qr(x)$rank)
}

evaluate_estimates <- function(train_true,
                               train_pred,
                               test_true,
                               test_pred,
                               beta = NULL, # for rank computation
                               gamma = NULL, # for rank computation
                               M = NULL, # for rank computation
                               model = "",
                               time = -1,
                               digits = 5){
  
  cbind(
  IMR::evaluate(train_pred, train_true) %>%
    dplyr::transmute(train_rmse = RMSE,
                     train_rrmse = Rel_RMSE,
                     train_spearman = Spearman_Rho),
  IMR::evaluate(test_pred, test_true) %>%
    dplyr::transmute(test_rmse = RMSE,
                     test_rrmse = Rel_RMSE,
                     test_spearman = Spearman_Rho)
  ) %>%
    round(digits) %>%
    mutate(time = time, model = model,
           rank_beta = ifelse(! is.null(beta) && is.matrix(beta), comp_rank(beta), NA),
           rank_gamma = ifelse(! is.null(gamma) && is.matrix(gamma), comp_rank(gamma), NA),
           sparse_beta = ifelse(! is.null(beta) && is.matrix(beta), mean(beta==0), NA),
           sparse_gamma = ifelse(! is.null(gamma) && is.matrix(gamma), mean(gamma==0), NA),
           rank_m = ifelse(! is.null(M) && is.matrix(M), comp_rank(M), NA))
}
