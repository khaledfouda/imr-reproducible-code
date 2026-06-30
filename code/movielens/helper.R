require(tidyverse)
require(magrittr)
library(IMR)
require(kableExtra)
require(RSSthemes)

#------------------------------------------------------------------------------
seed <- 2025
# configurations >>
source("./code/movielens/config_default.R")
if (file.exists("./code/movielens/config.R")) {
  source("./code/movielens/config.R")
}
# ===== loading and preparing the data ===============
# ===== out:   X, Z, Y, query, test.idx, test.truths, obs_mask
# =============================================================
load_movielens1m <- function() {
  data <- list()
  load("./data/movielens/raw/Movie_X.Rdata") # X
  load("./data/movielens/raw/Movie_Y.Rdata", verbose = T)
  X <- X[, 1:4] # keep only main-effects
  data$X <- X
  data$Y <- Y
  load("./data/movielens/raw/Movie_Q.Rdata")
  data$query <- query
  # ========================================================
  # prepare test set and X-QR
  data$test.idx <- cbind(data$query[, 1], data$query[, 2])
  data$test.truths <- data$query[, 3]
  data$Y <- IMR::as_incomplete(data$Y)
  data$obs_mask <- as.matrix((data$Y != 0) * 1)
  mean(data$obs_mask == 1)
  mean(data$obs_mask == 0)
  # ====================================================
  # prepare Z < genres >
  data$Z <- data.table::fread(
    file = "./data/movielens/raw/movies_Z.dat",
    sep = NULL,
    encoding = "Latin-1",
    header = FALSE
  ) %>%
    tidyr::separate(
      V1,
      into = c("movie_id", "title", "genres"),
      sep = "::"
    )

  genre_labels <- c(
    "Action", "Adventure", "Animation", "Children's", "Comedy", "Crime",
    "Documentary", "Drama", "Fantasy", "Film-Noir", "Horror", "Musical",
    "Mystery", "Romance", "Sci-Fi", "Thriller", "War", "Western"
  )
  data$Z %>%
    transmute(genre = as.vector(strsplit(data$Z$genres, "|", fixed = TRUE))) ->
  genres
  genres <- genres[[1L]]
  i <- rep(seq_along(genres), lengths(genres))
  j <- match(unlist(genres, use.names = FALSE), genre_labels)
  keep <- !is.na(j)
  genre_sparse <- Matrix::sparseMatrix(
    i = i[keep], j = j[keep], x = 1L,
    dims = c(length(genres), length(genre_labels)),
    dimnames = list(NULL, genre_labels)
  )
  genre_df <- as.data.frame(as.matrix(genre_sparse), check.names = FALSE)

  data$Z <- cbind(data$Z[, 1:2], genre_df)
  data$Z %<>% mutate(movie_id = as.numeric(movie_id))
  movies_no_genre <- (1:3952)[!(1:3952 %in% data$Z$movie_id)]
  extra.rows <- data.frame(movie_id = movies_no_genre, title = "")
  for (genre in genre_labels) {
    extra.rows[genre] <- 0
  }
  rbind(data$Z, extra.rows) %>%
    arrange(movie_id) %>%
    dplyr::select(-movie_id, -title) %>%
    as.matrix() ->
  data$Z
  return(data)
}
# =========
