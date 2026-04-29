source("./R/helper.R")
source("./R/other_models/SoftImpute_cv.R")
source("./R/other_models/MCCI.R")
source("./R/simulation/helper.R")

require(tidyverse)
require(IMR)

#==========================================
# setting 1)
#==========================================
dims <- c(400, 600, 800, 1000)
p = 4;
q = 0;
r = 10;
missing_pct = 0.8

all_res <- data.frame()

convergence <- IMR::imr_convergence(maxit=1000, thresh=1e-6)
grid <- IMR::imr_tune_grid(rank = c(2, 10, 1, 2), beta = 0, nuclear = c(0,40,40,2))
print(grid)

for(b in 1:500){
  seed = 2025 + b
  set.seed(seed)
  for(d in dims){
    
    n = m = d
    dat <-
      generate_simulated_data(n, m, r, p, q, missing_pct,
                              snr = 1,
                              shared = FALSE,
                              seed = seed)
    
    
    mdat <- IMR::imr_data(Y = dat$Y, X = dat$X, seed = seed, val_prop = 0.2);
    
    
    # fit the 3 models:
    
    fitsi <- simpute.cv(y_full = mdat$Y,
                        y_train = mdat$y_train,
                        y_valid = mdat$y_valid,
                        trace = FALSE,
                        print.best = FALSE,
                        tol = grid$nuclear$streaks,
                        n.lambda = grid$nuclear$length,
                        maxit = convergence$maxit,
                        thresh = convergence$thresh,
                        test_error = get_metric("rmse"),
                        seed = seed)
    
    fitimr <- IMR::imr_tune(mdat, grid, convergence=convergence, fast_nuclear = FALSE,
                            seed = seed, n_cores = 7, verbose = 0)
    
    
    fitmcci <- MCCI.cv(Y = dat$Y, X = dat$X, W = dat$mask, n_folds = 5,numCores = 9,
                       seed = seed,
                       test_error = get_metric("rmse"),
                       lambda_1_grid = c(0),#seq(0, 1, length = 10),
                       lambda_2_grid = seq(2.9, 0, length = 20),
                       alpha_grid = c(1),#seq(0.992, 1, length = 10),
                       n1n2_optimized = TRUE,
                       return_diagn = FALSE)
    #--- 
    # collect the results
    
    test_true <- dat$theta[dat$mask == 0]
    train_pred <- fitmcci$fit$estimates[dat$mask == 1]
    test_pred <- fitmcci$fit$estimates[dat$mask == 0]
    evaluate_estimates(mdat$Y@x, train_pred,
                       test_true, test_pred,
                       beta = fitmcci$fit$beta,
                       M = fitmcci$fit$M,
                       model = "MCCI") %>%
      mutate(beta_rrmse = evaluate(fitmcci$fit$beta, dat$beta, "rrmse"),
             M_rrmse = evaluate(fitmcci$fit$M, dat$M, "rrmse")) ->
      mcci_out
    

    test_mask <- dat$mask == 0
    test_mask <- as_incomplete(test_mask)
    reconst <- reconstruct(fitimr$fit, mdat, FALSE)
    
    train_pred <- reconstruct_partial(fitimr$fit, mdat, mdat$Y@i, mdat$Y@p)
    test_pred <- reconstruct_partial(fitimr$fit, mdat, test_mask@i, test_mask@p)
    evaluate_estimates(mdat$Y@x, train_pred,
                       test_true, test_pred,
                       beta = reconst$beta,
                       M = reconst$M,
                       model = "IMR") %>%
      mutate(beta_rrmse = evaluate(reconst$beta, dat$beta, "rrmse"),
             M_rrmse = evaluate(reconst$M, dat$M, "rrmse")) ->
      imr_out
    
    e <- fitsi$fit
    estimates <- e$u %*% (t(e$v) * e$d)
    train_pred <- estimates[dat$mask == 1]
    test_pred <- estimates[dat$mask == 0]
    evaluate_estimates(mdat$Y@x, train_pred,
                       test_true, test_pred,
                       M = estimates,
                       model = "SI") %>%
      mutate(beta_rrmse = NA, M_rrmse = NA) ->
      si_out
    
    rbind(mcci_out, imr_out, si_out) %>%
      mutate(dim = d,
             rank = rank_m +  if_else(is.na(rank_beta), 0, rank_beta)) ->
      result
    all_res <- rbind(all_res, result)
  }
  print(b)
  all_res %>%
    group_by(model, dim) %>%
    summarize_all(c(m=mean,s=sd)) %>%
    as.data.frame() %>%
    select(model, dim, test_rrmse_m, test_rrmse_s, rank_m) %>%
    arrange(dim, test_rrmse_m) %>%
    print()
  
  rw_a_file("results_scenario_1.rds",
            data = all_res,
            file_override = TRUE,
            directory = "./data/Simulation/",
            type = "write"
            )
  }
#============================================================
# Setting 2 - part 1 (to generate Figure 1) 
# Models: SImpute and IMR
#============================================================
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

convergence <- IMR::imr_convergence(maxit=1000, thresh=1e-6)
grid <- IMR::imr_tune_grid(rank = c(2, 10, 1, 2), beta=c(0), gamma=c(0), nuclear=c(0,120,60,2));
print(grid)

n = m = 1000
p = 5;
q = 5;
r = 5;
missing_pct = seq(.7, .98, .05)
all_res <- res <- data.frame()
# b = 1; pct=1
for(b in 1:500){
  seed = 2025 + b
  start1 = Sys.time()
  set.seed(seed)
  dat <-
    generate_simulated_data(n, m, r, p, q, .7,
                            snr = 1,
                            seed = seed)
  for(pct in 1:length(missing_pct)){
    start2 = Sys.time()
    if(pct > 1)
      dat <- increase_sparsity(dat, .05)
    
    mdat <- IMR::imr_data(Y = dat$Y, X = dat$X, Z = dat$Z,  seed = seed, val_prop = 0.2);
    
    # fit the 3 models:
    
    fitsi <- simpute.cv(y_full = mdat$Y,
                        y_train = mdat$y_train,
                        y_valid = mdat$y_valid,
                        trace = FALSE,
                        print.best = FALSE,
                        tol = grid$nuclear$streaks,
                        maxit = convergence$maxit,
                        thresh = convergence$thresh,
                        n.lambda = grid$nuclear$length,
                        test_error = get_metric("rmse"),
                        seed = seed)
    
    fitimr <- IMR::imr_tune(mdat, grid, convergence=convergence, fast_nuclear = FALSE,
                            seed = seed, n_cores = 7, verbose = 0)
    
    
    #--- 
    # collect the results
    
    test_true <- dat$theta[dat$mask == 0]
    
    test_mask <- dat$mask == 0
    test_mask <- as_incomplete(test_mask)
    reconst <- reconstruct(fitimr$fit, mdat, FALSE)
    
    train_pred <- reconstruct_partial(fitimr$fit, mdat, mdat$Y@i, mdat$Y@p)
    test_pred <- reconstruct_partial(fitimr$fit, mdat, test_mask@i, test_mask@p)
    evaluate_estimates(mdat$Y@x, train_pred,
                       test_true, test_pred,
                       beta = reconst$beta,
                       gamma = reconst$gamma,
                       M = reconst$M,
                       model = "IMR") %>%
      mutate(lambda_m = fitimr$params$lambda_m,
             rank_in = fitimr$params$rank_in,) %>%
      mutate(beta_rrmse = evaluate(reconst$beta, dat$beta, "rrmse"),
             gamma_rrmse = evaluate(reconst$gamma, dat$gamma, "rrmse"),
             M_rrmse = evaluate(reconst$M, dat$M, "rrmse"),
             theta_rrmse = evaluate(reconst$estimates, dat$theta, "rrmse")) ->
      imr_out
    
    e <- fitsi$fit
    estimates <- e$u %*% (t(e$v) * e$d)
    train_pred <- estimates[dat$mask == 1]
    test_pred <- estimates[dat$mask == 0]
    evaluate_estimates(mdat$Y@x, train_pred,
                       test_true, test_pred,
                       M = estimates,
                       model = "SI") %>%
      mutate(lambda_m = fitsi$lambda,
             rank_in = fitsi$rank.max, beta_rrmse = NA, gamma_rrmse = NA,
             M_rrmse = NA,
             theta_rrmse = evaluate(estimates, dat$theta, "rrmse")) ->
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
    
  }
  all_res %>%
    group_by(model, missing_pct) %>%
    summarize_all(c(m=mean,s=sd)) %>%
    as.data.frame() %>%
    select(model, missing_pct, test_rrmse_m, test_rrmse_s, rank_m, lambda_m_m, rank_in_m) %>%
    arrange(missing_pct, test_rrmse_m) %>%
    print()
  
  rw_a_file("results_scenario_2_part_1.rds",
            data = all_res,
            file_override = TRUE,
            directory = "./data/Simulation/",
            type = "write"
  )
  print(paste(b, " in ", round(Sys.time() - start1)))
}
#============================================================
# Setting 2 - part 2 (to generate Figure 2) 
# Models: SImpute and IMR
#============================================================
results_scenario2_part_1 <- rw_a_file("results_scenario_2_part_1_backup.rds",
                                      directory = "./data/Simulation/",
                                      type = "read"
)

# extract hyperparameters: choose the median of overall
results_scenario2_part_1 %>%
  filter(model == "IMR") %>%
  dplyr::select(lambda_laplace, rank_m, miss_pct) %>%
  mutate(miss_pct = round(miss_pct, 2)) %>%
  summarise_all(median) %>%
  mutate(rank_m = round(rank_m)) %>%
  ungroup() -> best_hparams

results_scenario2_part_1 %>%
  filter(model == "SI") %>%
  dplyr::select(rank, miss_pct, lambda_laplace) %>%
  mutate(miss_pct = round(miss_pct, 2)) %>%
  summarise_all(median) %>%
  round() -> simpute_rank
#----------------------------------------------------------------

n = m = 1000
p = 5;
q = 5;
r = 5;
missing_pct = seq(.7, .98, .05)
convergence <- IMR::imr_convergence(maxit=1000, thresh=1e-6)
all_res <- res <- data.frame()

for(b in 1:500){
  seed = 2025 + b
  start1 = Sys.time()
  set.seed(seed)
  dat <-
    generate_simulated_data(n, m, r, p, q, .7,
                            snr = 1,
                            seed = seed)
  for(pct in 1:length(missing_pct)){
    start2 = Sys.time()
    if(pct > 1)
      dat <- increase_sparsity(dat, .05)
    
    mdat <- IMR::imr_data(Y = dat$Y, X = dat$X, Z = dat$Z,  seed = seed, val_prop = 0.2)
    
    start = Sys.time()
    fitsi <- softImpute::softImpute(dat$Y,
                                    rank.max = simpute_rank$rank,
                                    lambda = simpute_rank$lambda_laplace,
                                    thresh = convergence$thresh,
                                    maxit = convergence$maxit,
                                    trace.it = FALSE,final.svd = TRUE, type = "als")
    time.si = as.numeric(Sys.time() -  start, units = "secs")
    
    start = Sys.time()
    fitimr <- IMR::imr_fit(mdat, rank = best_hparams$rank_m,
                           lambda_m = best_hparams$lambda_laplace,
                           convergence=convergence)
    time.imr = as.numeric(Sys.time() -  start, units = "secs")
    
    
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
  }
  print(paste(b, " in ", round(Sys.time() - start1)))
  
  all_res %>%
    mutate(missing_pct = round(missing_pct,2)) %>%
    group_by(model, missing_pct) %>%
    summarize_all(c(m=mean,s=sd)) %>%
    as.data.frame() %>%
    select(model, missing_pct, test_rrmse_m, test_rrmse_s, rank_m) %>%
    arrange(missing_pct, test_rrmse_m) %>%
    print()
  
  rw_a_file("results_scenario_2_part_2.rds",
            data = all_res,
            file_override = TRUE,
            directory = "./data/Simulation/",
            type = "write"
  )
}
