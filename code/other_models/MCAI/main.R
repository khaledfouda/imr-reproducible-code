# source of files: 
# Supplementary material of
# Statistical Inference For Noisy Matrix Completion Incorporating Auxiliary Information. Taylor & Francis. 
# Ma, Shujie; Niu, Po-Yao; Zhang, Yichong; Zhu, Yinchu.


source('./code/other_models/MCAI/utilities.r')
source('./code/other_models/MCAI/MC_solver.r')
source('./code/other_models/MCAI/MC_LS_class.r') 



MCAI.fit <- function(Y,
                     X,
                     max_iter = 30L,
                     tol=1e-6,
                     C_h=0.2,
                     delta_h=0.1,
                     r_bar = 8,
                     save_to_file=FALSE,
                     file_location = NULL){
  
  start_time <- Sys.time()
  M <- MC_alt_LS_rank_solver(Y=Y, X=X, intercept_val=1,
                             max_iter = max_iter, .tol = tol, missing_model = 'logistic')
  M$.make_obs()
  # estimate missing model
  M$estimate_pi()
  #R <- summary(M$pi_est$fit)
  #Tab <- R$coefficients
  #rownames(Tab) <- c("(Intercept)", colnames(X))
  
  # Include only main effects for model fitting (best performing based on the article)
  M$set_fitting_cov(index = 1:(ncol(X)+1) )
  print(M$summary())
  # Rank estimation
  M$rank_estimation(C_h=C_h, delta_h=delta_h, penalty_names = c('h'), r_bar=r_bar)
  rhat <- M$rank_est$est['h'] # get the rank estimation result
  # Fitting with the estimated rank
  M$fitting(with_rs = rhat)
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  out <- list(fit = M,
              time_secs = total_time)
  
  if(save_to_file){
    if(is.null(file_location)){
      message("file location must be provided")
    }else{
      message(paste("Saving file at :", file_location))
      saveRDS(out, file_location)
    }
  }
  return(out)
}
