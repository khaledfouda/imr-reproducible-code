source("./R/helper.R")
source("./R/Bixi/helper.R")
library(IMR)

seed <- 4000
convergence <- imr_convergence(maxit = 2000,
                               thresh = 1e-7)
#-------------------------------------------------------
# fit all 50 * 5 data files.
train_seq <- seq(55, 75, 5)
num_splits <- 50

split_id <- 1
train_size <- 55
total_results <- data.frame()

split_id = 1; train_size = train_seq[1]; # used for diagnostics

for(split_id in 1:num_splits){
  for(train_size in train_seq){
    

bixi <- bixi_load_split(split_id, train_size)

kernels <- bixi_load_kernels(split_id, train_size)

sim_rows <- imr_similarity(kernels$temporal, jitter = 0.5, invert = TRUE); sim_rows
sim_cols <- imr_similarity(kernels$spatial, jitter = 0.5, invert = TRUE); sim_cols

data <- imr_data(bixi$Y, similarity_rows = sim_rows, similarity_cols = sim_cols,
                 val_prop = 0.2, seed = seed)
# 1. fit the model with similarity and intercepts
data <- update(data, row_intercept = TRUE, col_intercept  = TRUE)

start_time <- Sys.time()
fitsim <- imr_fit(
  data = data,
  rank = 12,
  lambda_m = 1.69,
  convergence = convergence
)
total_time_secs <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
output <- reconstruct(fitsim, data,FALSE)

train_pred <- reconstruct_partial(fitsim, data, data$Y@i, data$Y@p)
test_pred <- reconstruct_partial(fitsim, data, bixi$test@i, bixi$test@p)

total_results <- rbind(total_results,
    evaluate_estimates(data$Y@x, train_pred,
                       bixi$test@x, test_pred,
                       time = total_time_secs,
                       model = "IMR-S") %>%
      mutate(rank_m = sum(fitsim$coefficients$d > 0),
             split_id = split_id,
             train_size = train_size))

# 2. fit the model with intercepts only
data <- update(data, row_similarity = FALSE, col_similarity = FALSE)

start_time <- Sys.time()
fitn <- imr_fit(
  data = data,
  rank = 11,
  lambda_m = 0.775,
  convergence = convergence
)
total_time_secs <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
output <- reconstruct(fitn, data,FALSE)

train_pred <- reconstruct_partial(fitn, data, data$Y@i, data$Y@p)
test_pred <- reconstruct_partial(fitn, data, bixi$test@i, bixi$test@p)

total_results <- rbind(total_results,
                      evaluate_estimates(data$Y@x, train_pred,
                                         bixi$test@x, test_pred,
                                         time = total_time_secs,
                                         model = "IMR-N") %>%
                        mutate(rank_m = sum(fitn$coefficients$d > 0),
                               split_id = split_id,
                               train_size = train_size))


  }
}

# save to disk >> 
rw_a_file(
  "results_imr.rds",
  data = total_results,
  file_override = TRUE,
  create_folder = TRUE,
  directory = "./data/Bixi/",
  type = "write"
)


# look at the results
total_results %>%
  group_by(model, train_size) %>%
  summarize_all(mean) %>%
  group_by(model) %>%
  mutate_all(round, 4) %>%
  ungroup() %>%
  as.data.frame() %>%
  arrange(train_size, test_rmse) %>% 
  print()
