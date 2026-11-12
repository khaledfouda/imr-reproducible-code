export_movielens_for_glocalk <- function(
    y_rdata_path = "./data/movielens/raw/Movie_Y.Rdata",
    q_rdata_path = "./data/movielens/raw/Movie_Q.Rdata",
    y_out_path = "./data/movielens/raw/Movie_Y_Glocalk.dat",
    test_out_path = "./data/movielens/raw/Movie_test_Glocalk.dat"
) {
  require(Matrix)
  message("Loading Rdata files...")
  load(y_rdata_path) # Loads 'Y'
  load(q_rdata_path) # Loads 'query'
  
  message("Formatting Y matrix...")
  obs_ind <- which(Y != 0, arr.ind = TRUE)
  py.Y <- data.frame(userID = obs_ind[, 1], movieID = obs_ind[, 2], rating = Y[obs_ind])
  
  message("Formatting query matrix...")
  query <- as.data.frame(query)
  colnames(query) <- c("userID", "movieID", "rating")
  
  message("Writing Y to ", y_out_path)
  write.table(py.Y,
              y_out_path,
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
              
  message("Writing test to ", test_out_path)
  write.table(query,
              test_out_path,
              sep       = "::",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
              
  message("Done! Files are ready for Glocal-K.")
}

export_movielens_for_glocalk()
