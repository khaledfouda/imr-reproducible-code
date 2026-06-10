source("./R/helper.R")
source("./R/other_models/SoftImpute_cv.R")
source("./R/other_models/MCCI.R")
source("./R/simulation/helper.R")

require(tidyverse)
require(IMR)

#==========================================
# Initial settings copied from 1_simulations.R
#==========================================
dims <- c(400, 600, 800, 1000)
p = 4;
q = 0;
r = 10;
missing_pct = 0.8
increase_sparsity <- function(dat, step=0.05){
  current_sparsity <- mean(dat$mask == 0)
  target_sparsity <- step + current_sparsity
  print(paste("Target sparsity is ", target_sparsity))
  stopifnot(target_sparsity < 1)
  extra_nonzero_frac <- step / (1 - current_sparsity)
  nonzero_idx <- which(dat$mask == 1)
  to_zero_ind <- sample(nonzero_idx, extra_nonzero_frac*length(nonzero_idx),replace = F)
  
  dat$Y[to_zero_ind] <- 0
  #dat$Y %<>% IMR::as_incomplete()
  #-- we now recreate the train/test splits
  dat$mask <- as.matrix(dat$Y != 0)
  dat$sparsity <- target_sparsity
  return(dat)
}


results_scenario2_part_1 <- rw_a_file("results_scenario_2_part_1.rds",
                                      directory = "./data/Simulation/",
                                      type = "read"
)

# extract hyperparameters: choose the median of overall
results_scenario2_part_1 %>%
  filter(model == "IMR") %>%
  dplyr::select(lambda_m, rank_m, missing_pct) %>%
  mutate(missing_pct = round(missing_pct, 2)) %>%
  summarise_all(median) %>%
  mutate(rank_m = round(rank_m)) %>%
  ungroup() -> best_hparams

results_scenario2_part_1 %>%
  filter(model == "SI") %>%
  dplyr::select(rank, missing_pct, lambda_m) %>%
  mutate(missing_pct = round(missing_pct, 2)) %>%
  summarise_all(median) %>%
  round() -> simpute_rank
#----------------------------------------------------------------

n = m = 1000
p = 5;
q = 5;
r = 5;
missing_pct = seq(.7, .98, .05)
convergence <- IMR::imr_convergence(maxit=1000, thresh=1e-6, trace=TRUE)
all_res <- res <- data.frame()

# run a single iteration
b=1; pct = 1

seed = 2025 + b
set.seed(seed)
#-- prepare the data
dat <-
  generate_simulated_data(n, m, r, p, q, .7,
                          snr = 1,
                          seed = seed)
if(pct > 1)
  dat <- increase_sparsity(dat, .05)

mdat <- IMR::imr_data(Y = dat$Y, X = dat$X, Z = dat$Z,  seed = seed, val_prop = 0.2)

# fit SI
start = Sys.time()
fitsi <- softImpute::softImpute(dat$Y,
                                rank.max = simpute_rank$rank,
                                lambda = simpute_rank$lambda_m,
                                thresh = convergence$thresh,
                                maxit = convergence$maxit,
                                trace.it = TRUE,final.svd = TRUE, type = "als")
time.si = as.numeric(Sys.time() -  start, units = "secs")


# fit IMR
start = Sys.time()
fitimr <- IMR::imr_fit(mdat, rank = best_hparams$rank_m,
                       lambda_m = best_hparams$lambda_m,
                       convergence=convergence)
time.imr = as.numeric(Sys.time() -  start, units = "secs")

f





#--- 
# collect the results

test_true <- dat$theta[dat$mask == 0]

test_mask <- dat$mask == 0
test_mask <- as_incomplete(test_mask)
reconst <- reconstruct(fitimr, mdat, FALSE)

train_pred <- reconstruct_partial(fitimr, mdat, mdat$Y@i, mdat$Y@p)
test_pred <- reconstruct_partial(fitimr, mdat, test_mask@i, test_mask@p)
evaluate_estimates(mdat$Y@x, train_pred,
                   test_true, test_pred,
                   beta = reconst$beta,
                   gamma = reconst$gamma,
                   M = reconst$M,
                   time = time.imr,
                   model = "IMR")  ->
  imr_out

e <- fitsi
estimates <- e$u %*% (t(e$v) * e$d)
train_pred <- estimates[dat$mask == 1]
test_pred <- estimates[dat$mask == 0]
evaluate_estimates(mdat$Y@x, train_pred,
                   test_true, test_pred,
                   M = estimates,
                   time = time.si,
                   model = "SI")  ->
  si_out

rbind(imr_out, si_out) %>%
  mutate(dim = n,
         p = p,
         r = r,
         q = q,
         missing_pct = mean(dat$mask == 0),
         rank = rank_m +  if_else(is.na(rank_beta), 0, rank_beta) +
           if_else(is.na(rank_gamma), 0, rank_gamma)) ->
  result
all_res <- rbind(all_res, result)
