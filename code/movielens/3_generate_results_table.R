# =====================================================================
# MovieLens 1M — Analysis & Visualisation
# =====================================================================
# Produces: results comparison table + boxplot of fitted ratings
# Requires: model fits from MovieLens_runs.R
# =====================================================================

source("./code/helper.R")
source("./code/movielens/helper.R")

# ---  Load data --------------------------------------------------------
data <- load_movielens1m()
seed <- 2025
model_data <- imr_data(data$Y, data$X, data$Z, seed = seed, val_prop = 0)
print(model_data)

# --- Load fitted models -----------------------------------------------
imr_i <- readRDS("./output/movielens/model_fits/IMR_I_fit.rds")
imr_ixz <- readRDS("./output/movielens/model_fits/IMR_IXZ_fit.rds")
simpute <- readRDS("./output/movielens/model_fits/simpute_fit.rds")
mcai <- readRDS("./output/movielens/model_fits/mcai_fit.rds")

# ---  Reconstruct estimates --------------------------------------------

#  IMR models
imr_i$rank_M <- imr_i$meta$rank
imr_ixz$rank_M <- imr_ixz$meta$rank
imr_i$time <- imr_i$time_secs / 60
imr_ixz$time <- imr_ixz$time_secs / 60

out <- list()
out[[1]] <- IMR::reconstruct(imr_i, model_data)
out[[2]] <- IMR::reconstruct(imr_ixz, model_data)

#  SoftImpute: reconstruct from SVD components (u, d, v)
simpute$time <- simpute$time_secs / 60
simpute$rank_M <- sum(simpute$d > 1e-5)
out[[3]] <- list(
  estimates = simpute$u %*% (simpute$d * t(simpute$v)),
  beta = NULL,
  gamma = NULL
)

#  MCAI (Ma et al. 2025)
fmcai <- mcai$fit$fit[[1]]
mcai$rank_M <- fmcai$rank
mcai$time <- mcai$time_secs / 60
fmcai$M <- fmcai$L %*% t(fmcai$R)
fmcai$xbeta <- cbind(1, data$X) %*% t(fmcai$B)
fmcai$estimates <- fmcai$M + fmcai$xbeta
fmcai$beta <- t(fmcai$B[, -1])
out[[4]] <- fmcai
mcai$fit <- fmcai

# ---  Evaluate all models ----------------------------------------------

fits <- list(imr_i, imr_ixz, simpute, mcai)
models <- c("IMR-I", "IMR-IXZ", "SoftImpute", "MCAI")

res <- data.frame()
for (i in seq_along(models)) {
  message(models[[i]])
  res <- rbind(
    res,
    evaluate_estimates(
      model_data$Y[model_data$Y != 0],
      as_incomplete(out[[i]]$estimates * data$obs_mask)@x,
      data$test.truths,
      out[[i]]$estimates[data$test.idx],
      beta = out[[i]]$beta,
      gamma = out[[i]]$gamma,
      time = fits[[i]]$time,
      model = models[[i]]
    ) %>% mutate(rank_m = fits[[i]]$rank_M)
  )
}


res <- res %>%
  mutate(
    rank_beta  = rank_beta + 1, # extra rank for row intercepts
    rank_gamma = rank_gamma + 1, # extra rank for column intercepts
    rank_beta  = if_else(model == "IMR-I", 1, rank_beta),
    rank_gamma = if_else(model == "IMR-I", 1, rank_gamma)
  ) |>
  rbind(
    read.csv("./output/movielens/model_fits/glocalk_results.csv")
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(rank_total = rank_m +
    ifelse(is.na(rank_beta), 0, rank_beta) +
    ifelse(is.na(rank_gamma), 0, rank_gamma)) |>
  dplyr::arrange(test_rmse) |>
  transmute(
    model, time, rank_beta, rank_gamma, rank_total,
    error.test = test_rmse,
    rank_M = rank_m,
    corr.test = test_spearman,
    error.train = train_rmse,
    sparsity_beta = sparse_beta,
    sparsity_gamma = sparse_gamma
  ) -> res_df

# ---  Results table  --------------------------------------------

best_min_cols <- c("time", "error.test", "error.train", "rank_M")
best_max_cols <- c("corr.test")

best_idx_min <- lapply(best_min_cols, function(v) {
  which.min(replace(res_df[[v]], is.na(res_df[[v]]), Inf))
}) |> rlang::set_names(best_min_cols)

best_idx_max <- lapply(best_max_cols, function(v) {
  which.max(replace(res_df[[v]], is.na(res_df[[v]]) | is.nan(res_df[[v]]), -Inf))
}) |> rlang::set_names(best_max_cols)

fmt_num <- function(x, digits) {
  ifelse(is.na(x) | is.nan(x), "\u2014", sprintf(paste0("%.", digits, "f"), x))
}

disp <- res_df |>
  dplyr::mutate(
    time           = fmt_num(time, 2),
    error.test     = fmt_num(error.test, 3),
    corr.test      = fmt_num(corr.test, 3),
    error.train    = fmt_num(error.train, 3),
    sparsity_beta  = fmt_num(sparsity_beta, 3),
    sparsity_gamma = fmt_num(sparsity_gamma, 3),
    rank_M         = ifelse(is.na(rank_M), "\u2014", as.character(rank_M)),
    rank_beta      = ifelse(is.na(rank_beta), "\u2014", as.character(rank_beta)),
    rank_gamma     = ifelse(is.na(rank_gamma), "\u2014", as.character(rank_gamma)),
    rank_total     = as.character(rank_total)
  ) |>
  dplyr::select(-rank_total)

for (v in names(best_idx_min)) {
  idx <- best_idx_min[[v]]
  disp[[v]] <- kableExtra::cell_spec(disp[[v]], "latex", bold = seq_len(nrow(disp)) == idx)
}

for (v in names(best_idx_max)) {
  idx <- best_idx_max[[v]]
  disp[[v]] <- kableExtra::cell_spec(disp[[v]], "latex", bold = seq_len(nrow(disp)) == idx)
}


kbl(
  disp,
  format = "simple",
  booktabs = TRUE,
  linesep = "",
  escape = FALSE,
  caption = paste0(
    "Performance comparison on the MovieLens 1M dataset. ",
    "Best values per column are bolded and IMR models are shaded."
  )
) |>
  add_header_above(c(" " = 2, "Test" = 2, "Train" = 1, " " = 5)) |>
  add_header_above(c(" " = 2, "Performance" = 3, "Rank estimation" = 3, "Sparsity" = 2)) |>
  kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 8) |>
  row_spec(which(grepl("^IMR", disp$model)), bold = FALSE, background = "#f7f7f7") |>
  column_spec(1, width = "3.2cm")

# ===  Boxplot of fitted ratings ========================================

# Movie-genre mapping for selected genres
genres <- c(
  "Documentary", "Musical", "Drama", "Fantasy", "Children's",
  "War", "Action", "Sci-Fi", "Horror", "Animation"
)

z_sel <- data$Z[, genres, drop = FALSE]
arr <- which(z_sel != 0, arr.ind = TRUE)

map_z <- data.frame(
  genre = colnames(z_sel)[arr[, 2]],
  movie = as.character(arr[, 1]),
  stringsAsFactors = FALSE
) |> dplyr::mutate(movie = as.numeric(movie))

# User demographic groups
map_x <- as.data.frame(data$X) |>
  dplyr::mutate(user = seq_len(nrow(data$X))) |>
  dplyr::mutate(
    group = dplyr::case_when(
      .G == 1 & .A1 == 1 ~ "Male (25-34)",
      .G == 0 & .A1 == 1 ~ "Female (25-34)",
      .G == 1 & .A2 == 1 ~ "Male (35-49)",
      .G == 0 & .A2 == 1 ~ "Female (35-49)",
      .G == 1 & .A3 == 1 ~ "Male (50+)",
      .G == 0 & .A3 == 1 ~ "Female (50+)",
      .G == 1 ~ "Male (0-24)",
      .G == 0 ~ "Female (0-24)"
    )
  ) |>
  dplyr::select(user, group)

# Keep only movies belonging to a single genre
map_z <- map_z |>
  dplyr::group_by(movie) |>
  dplyr::mutate(n = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::filter(n == 1) |>
  dplyr::select(-n)

# Extract fitted estimates from IMR-IXZ
movies <- unique(arr[, 1])
sub_yh <- out[[2]]$estimates[, movies]
colnames(sub_yh) <- movies

# Long-format user-movie-estimate table
ume <- as.data.frame(sub_yh) |>
  tibble::rownames_to_column("row") |>
  tidyr::pivot_longer(-row, names_to = "col", values_to = "value") |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
  dplyr::rename(movie = col, user = row, estimate = value) |>
  dplyr::inner_join(map_z, by = "movie", relationship = "many-to-one") |>
  dplyr::inner_join(map_x, by = "user", relationship = "many-to-one")

q7 <- function(x, p) {
  stats::quantile(x, probs = p, na.rm = TRUE, type = 7, names = FALSE)
}

# Compute per-movie IQR
pooled_iqr_tbl <- ume |>
  dplyr::summarise(
    q1 = q7(estimate, 0.25),
    q3 = q7(estimate, 0.75),
    pooled_iqr = q3 - q1,
    .by = c(movie, genre)
  ) |>
  dplyr::select(movie, genre, pooled_iqr)

# Compute per-movie gap between group medians
median_gap_tbl <- ume |>
  dplyr::summarise(
    med = median(estimate, na.rm = TRUE),
    .by = c(movie, genre, group)
  ) |>
  dplyr::summarise(
    s_diff = diff(range(med)),
    .by = c(movie, genre)
  )

# Movie titles
titles <- data.table::fread(
  file = "./data/movielens/raw/movies_Z.dat",
  sep = NULL,
  encoding = "Latin-1",
  header = FALSE
) |>
  tidyr::separate(V1, into = c("movie", "title", "genres"), sep = "::") |>
  as.data.frame() |>
  dplyr::mutate(movie = as.numeric(movie)) |>
  dplyr::select(-genres) |>
  dplyr::filter(movie %in% map_z$movie)

# Select 3 children's movies with highest heterogeneity score
movies_picked <- c(586, 1367, 1592)
# "Home Alone (1990)", "101 Dalmatians (1996)", "Air Bud (1997)"

selected_movies <- pooled_iqr_tbl |>
  dplyr::inner_join(median_gap_tbl, by = c("movie", "genre")) |>
  dplyr::mutate(
    .by = genre,
    z_iqr = as.numeric(scale(pooled_iqr)),
    z_diff = as.numeric(scale(s_diff)),
    s = 0.5 * z_diff + 0.5 * z_iqr
  ) |>
  dplyr::filter(movie %in% movies_picked, genre == "Children's") |>
  dplyr::arrange(genre, dplyr::desc(s))

selected_movies <- titles |>
  dplyr::filter(movie %in% selected_movies$movie) |>
  dplyr::inner_join(selected_movies, by = "movie", relationship = "one-to-one") |>
  dplyr::select(-z_iqr, -z_diff)

ume_sel <- ume |>
  dplyr::semi_join(selected_movies, by = c("movie", "genre")) |>
  dplyr::inner_join(titles, by = "movie", relationship = "many-to-one")

group_lv <- c(
  "Female (0-24)", "Female (25-34)", "Female (35-49)", "Female (50+)",
  "Male (0-24)",   "Male (25-34)",   "Male (35-49)",   "Male (50+)"
)
age_lv <- c("0-24", "25-34", "35-49", "50+")

ume_sel <- ume_sel |>
  dplyr::mutate(
    group  = factor(group, levels = group_lv),
    age    = stringr::str_extract(group, "(?<=\\().+(?=\\))") |> factor(levels = age_lv),
    gender = ifelse(stringr::str_detect(group, "^Female.*"), "Female", "Male") |> factor()
  )

# Order movies within genre by selection score
movie_order <- selected_movies |>
  dplyr::arrange(genre, dplyr::desc(s)) |>
  dplyr::distinct(genre, movie) |>
  dplyr::group_by(genre) |>
  dplyr::mutate(order = dplyr::row_number()) |>
  dplyr::ungroup()

ume_sel <- ume_sel |>
  dplyr::left_join(movie_order, by = c("movie", "genre"))

ume_sel$movie <- factor(
  ume_sel$movie,
  levels = unique(ume_sel$movie[order(ume_sel$genre, ume_sel$order)])
)

# Build the plot
library(ggh4x)

g <- ume_sel %>%
  rename(`Age Group` = age) %>%
  ggplot(aes(x = `Age Group`, y = estimate, fill = `Age Group`)) +
  geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 0.7, alpha = 0.9) +
  ggh4x::facet_nested(
    cols = vars(title, gender),
    scales = "fixed",
    strip = ggh4x::strip_nested(
      text_x = ggh4x::elem_list_text(face = c("bold"))
    )
  ) +
  scale_fill_rss_d(palette = "signif_seq", -1) +
  scale_y_continuous(
    "Fitted Rating",
    limits = c(0.5, 5),
    breaks = seq(0, 5, 1)
  ) +
  theme_significance() +
  theme(
    legend.position = "none",
    legend.justification = "left",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey98")
  )

print(g)

# Save plot
ggsave(
  filename = "./output/movielens/plot_full_model.pdf",
  plot = g,
  device = cairo_pdf,
  width = 12,
  height = 6,
  units = "in"
)
