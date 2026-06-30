source("./code/helper.R")
source("./code/movielens/helper.R")
source("./code/other_models/SoftImpute_cv.R")
source("./code/other_models/MCAI/main.R")

# if true, IMR will run cross-validation.
# load the data
data <- load_movielens1m()
model_data <- imr_data(data$Y, data$X, data$Z, seed = seed, val_prop = 0.2)
print(model_data)

# convergence for tuning and convergence2 for fitting.
convergence <- imr_convergence(maxit = 1000, thresh = 1e-5)
convergence2 <- imr_convergence(maxit = 1000, thresh = 1e-7)

# grid of hyperparameters
grid <- imr_tune_grid(
  beta = c(0, 0.4, 60),
  gamma = c(0, 0.4, 60),
  nuclear = c(0, 20, 80, 3),
  rank = c(5, 20, 1, 3)
)


# 1. fit imr with covariates

model_data <- update(model_data,
  row_intercept = TRUE,
  col_intercept = TRUE
)
model_data

# no need to set grid limits since we decided on the max values
# grid <- imr_set_grid_limits(model_data, grid,default_rank = 10,
#                             bisection_iter = 5,
#                             convergence=convergence, verbose=2)

if (run_cross_validation) {
  imrxz_cv <- IMR::imr_tune(model_data, grid,
    final_fit = FALSE,
    fast_nuclear = FALSE,
    nuclear_log_scale = FALSE,
    convergence = convergence, n_cores = 7,
    seed = seed, verbose = 1
  )

  saveRDS(imrxz_cv, "./output/movielens/model_fits/IMR_IXZ_tune.rds")
} else {
  imrxz_cv <- readRDS("./output/movielens/model_fits/IMR_IXZ_tune.rds")
}

# we now fit using the optimal hyper-parameters.
start <- Sys.time()
imrxz_fit <- IMR::imr_fit(model_data,
  rank = imrxz_cv$params$rank,
  lambda_m = imrxz_cv$params$lambda_m,
  lambda_beta = imrxz_cv$params$lambda_beta,
  lambda_gamma = imrxz_cv$params$lambda_gamma,
  convergence = convergence2
)
time <- Sys.time() - start
imrxz_fit$time_secs <- as.numeric(time, "secs")
saveRDS(imrxz_fit, "./output/movielens/model_fits/IMR_IXZ_fit.rds")

print(imrxz_fit)
print(summary(imrxz_fit))
#------------------------------------------------------------------
# 2. fit imr without covariates

model_data <- update(model_data,
  row_covariates = FALSE,
  col_covariates = FALSE
)
model_data

# no need to set grid limits since we decided on the max values
# grid <- imr_set_grid_limits(model_data, grid,default_rank = 10,
#                             bisection_iter = 5,
#                             convergence=convergence, verbose=2)
if (run_cross_validation) {
  imri_cv <- IMR::imr_tune(model_data, grid,
    final_fit = FALSE,
    fast_nuclear = FALSE,
    nuclear_log_scale = FALSE,
    convergence = convergence, n_cores = 7,
    seed = seed, verbose = 1
  )
  saveRDS(imri_cv, "./output/movielens/model_fits/IMR_I_tune.rds")
} else {
  imri_cv <- readRDS("./output/movielens/model_fits/IMR_I_tune.rds")
}
# we now fit using the optimal hyper-parameters.
start <- Sys.time()
imri_fit <- IMR::imr_fit(model_data,
  rank = imri_cv$params$rank,
  lambda_m = imri_cv$params$lambda_m,
  lambda_beta = imri_cv$params$lambda_beta,
  lambda_gamma = imri_cv$params$lambda_gamma,
  convergence = convergence2
)
time <- Sys.time() - start
imri_fit$time_secs <- as.numeric(time, "secs")
saveRDS(imri_fit, "./output/movielens/model_fits/IMR_I_fit.rds")
print(imri_fit)
print(summary(imri_fit))

#------------------------------------------------------------------
# 3. fit Soft-Impute

# we begin by using cross-validation to select the hyperparameters
if (run_cross_validation) {
  start_si <- Sys.time()
  simpute_cv <- simpute.cv(
    y_full = model_data$Y,
    y_train = model_data$y_train,
    y_valid = model_data$y_valid,
    seed = seed,
    maxit = convergence$maxit,
    thresh = convergence$thresh,
    rank.limit = grid$rank$max,
    rank.step = grid$rank$step,
    rank.init = grid$rank$min,
    lambda_max = grid$nuclear$max,
    n.lambda = grid$nuclear$length,
    trace = TRUE
  )
  time_si <- Sys.time() - start_si
  simpute_cv$time_secs <- as.numeric(time_si, "secs")
  simpute_cv$time <- round(as.numeric(time_si, "mins"), 2)
  saveRDS(simpute_cv, "./output/movielens/model_fits/simpute_cv.rds")
} else {
  simpute_cv <- readRDS("./output/movielens/model_fits/simpute_cv.rds")
}
# we now fit using the chosen hyperparameters.
#  softimpute expects a normal matrix.
y_dense <- as.matrix(data$Y)
y_dense[y_dense == 0] <- NA
start <- Sys.time()
simpute_fit <- softImpute::softImpute(y_dense,
  rank.max = simpute_cv$rank.max,
  lambda = simpute_cv$lambda,
  thresh = convergence2$thresh,
  maxit = convergence2$maxit,
  trace.it = FALSE, final.svd = TRUE, type = "als"
)
simpute_fit$time_secs <- as.numeric(Sys.time() - start, units = "secs")
saveRDS(simpute_fit, "./output/movielens/model_fits/simpute_fit.rds")
#---------------------------------------------------------------------------------
# 4. we fit MCAI
load("./data/movielens/raw/Movie_Y.Rdata", verbose = T) # loads Y
M <- MCAI.fit(Y, data$X,
  save_to_file = T,
  file_location = "./output/movielens/model_fits/mcai_fit.rds"
)
#-------------------------------------------------------------------------------
