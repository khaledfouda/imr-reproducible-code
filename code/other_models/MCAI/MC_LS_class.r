#' @importClassesFrom Matrix Matrix
#' @importFrom methods setClassUnion setOldClass setRefClass
library(MASS)
library(Matrix)
library(methods)
library(magrittr)

# Update: 29 JAN 2024 ----- Niu, PoYao

setClassUnion("any_matrix", c("matrix", "Matrix"))

MC_model <- setRefClass('MC_model',
  # ===== Fields =====
  fields = list(
    # main data -----
    Y = 'any_matrix', 
    X = 'any_matrix', 
    
    # related to covariate X -----
    .colmeans_X = 'vector', 
    .colsds_X = 'vector', 
    intercept = 'numeric', 
    p = 'integer', 
    
    # related to response Y -----
    dim = 'vector', 
    .missing_type = 'character',
    missing_model = 'character',
    .obs = 'list',
    pi_est = 'list',
    obs_count = 'integer',
    
    # fitting parameters -----
    r = 'integer',
    fitting_cov='integer',
    pi_fitting_cov='integer',
    max_iter = 'integer',
    .tol = 'numeric',
    
    # fitting results -----
    ini_val = 'list',
    fit = 'ANY',
    rank_est = 'list',
    .solved = 'ANY',
    
    # user define -----
    other_info = 'list'
    
    # -----
  ),
  
  # ===== Methods =====
  methods = list(
    initialize = function(X=NULL, .missing_type = 'zero', missing_model = c('MCR', 'logistic')[1],
                          intercept_in_X = F, intercept_val = ifelse(intercept_in_X>0,1,0), 
                          max_iter = 30L, .tol = 1e-6, ...){
      
      .self$initFields(..., .missing_type = .missing_type, missing_model = missing_model, max_iter = max_iter, .tol = .tol)
      
      .self$dim <- dim(.self$Y)

      if(.self$.missing_type!='zero' && .self$.missing_type!='na') .self$.missing_type <- 'zero'
      
      .self$set_X(new_X = X, intercept_in_X = intercept_in_X, intercept_val = intercept_val)
    },
    
    # setting, modifying -----
    .make_obs = function(force = F){
      if (length(.obs)!=0 && !force){
        cat(sprintf('$.obs already exists ... (%s)\n', Sys.time()))
      } else {
        cat(sprintf('\nMaking obsevation index... (%s)\n', Sys.time()))
        if (.missing_type == 'na') tmp <- !is.na(Y)
        if (.missing_type == 'zero') tmp <- Y!=0

        .obs$byi <<- tmp %>% apply(1, which)
        .obs$byj <<- tmp %>% apply(2, which)
        .obs$all <<- which(tmp)
        cat(sprintf('$.obs done! ... (%s)\n', Sys.time()))
      }
      
      
      if (is.null(.obs$all)){
        if (.missing_type == 'na') tmp <- !is.na(Y)
        if (.missing_type == 'zero') tmp <- Y!=0
        .obs$all <<- which(tmp)
      }
      obs_count <<- length(.obs$all) %>% as.integer()
    },
    
    estimate_pi = function(forced = F){
      if (!is.null(pi_est$pi.hat) && !forced) {
        cat(sprintf('$pi_est$pi.hat already exists ... (%s)\n', Sys.time()))
        return(invisible(NULL))
      }
      .make_obs()
      cat(sprintf('\nEstimating pi with %s model ... (%s)\n', missing_model, Sys.time()))
      pi_est <<- pi.est(obs.index = .obs, miss.type = missing_model, X = X[,pi_fitting_cov,drop=F], detail = T)
      cat(sprintf('$pi_est done! ... (%s)\n', Sys.time()))
    },
    
    set_X = function(new_X=NULL, intercept_in_X = F, intercept_val = ifelse(intercept_in_X>0, 1, 0)){
      fitting_cov <<- integer(0)
      pi_fitting_cov <<- integer(0)
      if (is.null(new_X)) {
        .self$X <- matrix(0, nr = .self$dim[1], nc = 0)
      } else {
        k <- as.integer(intercept_in_X)
        if (k>0) .self$X <- new_X[,-k, drop=F] else .self$X <- new_X
        if (is.null(colnames(.self$X))) colnames(.self$X) <- paste0('V',1:ncol(.self$X))
      }
      intercept <<- 0
      .self$p <- dim(.self$X)[2]
      if (intercept_val!= intercept) .self$set_intercept(value = intercept_val, clear_results = F)
      .self$standardized_X <- F
      .self$.colmeans_X <- numeric()
      .self$.colsds_X <- numeric()
      if (p>0){
        W <- if(intercept==0) X else X[,-1, drop=F] 
        .colmeans_X <<- colMeans(W)
        .colsds_X <<- apply(W, MARGIN = 2,FUN = sd)
      }
      
      set_fitting_cov()
      set_pi_fitting_cov()
    },
    
    set_fitting_cov = function(index=NULL){
      if (ncol(X)==0){
        fitting_cov <<- integer(0)
        p <<- 0L
        clear_all_results(F)
        return(invisible(NULL))
      }
      if (is.null(index)) index <- 1:ncol(X)
      
      index <- sort(unique(as.integer(index[index<= ncol(X) & index>0])))
      if (length(index)==length(fitting_cov) && all(index==fitting_cov)) {
        cat('Same index. Nothing has changed.\n')
        return(invisible(NULL))
      }
      fitting_cov <<- index
      p <<- length(index)
      
      clear_all_results(F)
    },
    
    set_pi_fitting_cov = function(index=NULL){
      if (ncol(X)==0 || missing_model=='MCR'){
        pi_fitting_cov <<- integer(0)
        clear_pi_est()
        return(invisible(NULL))
      }
      if (is.null(index)) index <- 1:ncol(X)
      
      index <- sort(unique(as.integer(index[index<= ncol(X) & index>0])))
      if (length(index)==length(pi_fitting_cov) && all(index==pi_fitting_cov)) {
        cat('Same index. Nothing has changed.\n')
        return(invisible(NULL))
      }
      pi_fitting_cov <<- index

      clear_pi_est()
    },
    
    set_intercept  = function(value=1, clear_results = T){
      if (intercept == value) {
        cat('Same intercept value. Nothing has changed.\n')
        return(invisible(NULL))
      }
      if (intercept != 0 && value!=0) {
        X <<- cbind(value, X[,-1])
        if (1 %in% fitting_cov){
          if (length(fit)>0){
            for (k in 1:length(fit)){
              if (is.null(fit[[k]])) next
              fit[[k]]$B[,1] <<- fit[[k]]$B[,1]*intercept/value
            }
          }
          if (length(ini_val)>0){
            ini_val$B[,1] <<- ini_val$B[,1]*intercept/value
          }
        }
      } else if (intercept==0) {
        X <<- cbind(intercept=value, X)
        fitting_cov <<- fitting_cov + 1L
        pi_fitting_cov <<- pi_fitting_cov + 1L
      }else if (value==0) {
        X <<- X[,-1]
        if (1 %in% fitting_cov){
          clear_all_results(F)
        }
        fitting_cov <<- fitting_cov[fitting_cov>1] - 1L
        pi_fitting_cov <<- pi_fitting_cov[pi_fitting_cov>1] - 1L
      }
      intercept <<- value
    },
    
    # organizing -----
    clear_all_results = function(clear_pi_est=T){
      clear_initial()
      clear_rank_est()
      clear_fit()
      if (clear_pi_est) clear_pi_est()
    },
    
    clear_initial = function(){ 
      .self$ini_val <- list()
      cat('$ini_val: ')
      print(.self$ini_val)
    },
    
    clear_fit = function(clear_index = NULL){
      if (!is.null(clear_index)){
        fit[clear_index] <<- NULL
        keep_idx <- setdiff(1:length(r), clear_index)
        .solved <<- .solved[keep_idx]
        r <<- r[keep_idx]
      } else {
        fit <<- list()
        .solved <<- logical()
        r <<- integer()
        cat('$fit: ')
        print(.self$fit)
      }
    },
    
    clear_rank_est = function(){ 
      .self$rank_est <- list()
      cat('$rank_est: ')
      print(.self$rank_est)
    },
    
    clear_pi_est = function(){
      .self$pi_est <- list()
      cat('$pi_est: ')
      print(.self$pi_est)
    },
    
    rearrange_fit = function(ix = NULL){
      if (is.null(ix)) ix <- sort(r, index.return = T)$ix
      r <<- r[ix]
      fit <<- fit[ix]
      .solved <<- .solved[ix]
    }
    # -----
  )
)


# Basic class ----------
MC_alt_LS_solver <- setRefClass('MC_alt_LS_solver',
  contains = 'MC_model',
  
  # ===== Fields =====
  fields = list(
    # related to covariate X -----
    standardized_X = 'logical', 
    
    # related to response Y -----
    scaled_Y = 'logical',
    .multiplier_Y = 'numeric',
    .offset_Y = 'numeric'
    # -----
  ),
  
  # ===== Methods =====
  methods = list(
    
    initialize = function(max_iter = 3L, .tol = 1e-6, ...){
      callSuper(max_iter = max_iter, .tol = .tol, ...)
      .self$scaled_Y <- F
      .self$.solved <- rep(F, length(.self$r))
    },
    
    # pre processing -----
    scale_Y = function(stdize = T, colwise = T, off = NULL, mult = NULL){
      if (scaled_Y) stop('Already scaled: .$scaled_Y = TRUE')
      
      cat(sprintf('\nScaling Y matrix... (%s)\n', Sys.time()))
      m <- .self$dim[2]
      if (!is.null(off)) .offset_Y <<- if (colwise) rep(off, len=m) else off[1]
      if (!is.null(mult)) .multiplier_Y <<- if (colwise) rep(mult, len = m) else mult[1]
      
      
      find_offset <- function(v){
        if (diff(range(v))==0) {
          return(0)
        } else {
          return( ifelse(stdize, mean(v), min(v)) )
        }
      }
      
      find_multiplier <- function(v){
        drv <- diff(range(v))
        if (drv==0) {
          return(v[1])
        } else {
          return( ifelse(stdize, sd(v), drv) )
        }
      }
      
      .make_obs()
      
      if (!colwise){
        if (is.null(off)) .offset_Y <<- find_offset(Y[.obs$all])
        if (is.null(mult)) .multiplier_Y <<- find_multiplier(Y[.obs$all])
        Y[.obs$all] <<- (Y[.obs$all]-.offset_Y)/.multiplier_Y
      } else {
        if (is.null(off)) .offset_Y <<- numeric(m)
        if (is.null(mult)) .multiplier_Y <<- numeric(m)
        for (j in 1:m){
          i <- .obs$byj[[j]]
          if (length(i)==0) next
          if (is.null(off)) .offset_Y[j] <<- find_offset(Y[i, j])
          if (is.null(mult)) .multiplier_Y[j] <<- find_multiplier(Y[i, j])
          Y[i, j] <<- (Y[i, j]-.offset_Y[j])/.multiplier_Y[j]
        }
      }
      clear_all_results()
      scaled_Y <<- T
      
    },
    
    standardize_X = function(na.rm = T){
      if (p==0) stop('No X matrix\n')
      if (standardized_X) stop('Already standardized: .$standardized_X = TRUE\n')
      cat(sprintf('\nStandardizing X matrix... (%s)\n', Sys.time()))
      W <- if(intercept==0) X else X[,-1] 
      W <- t((t(W)-.colmeans_X)/.colsds_X)
      if (intercept != 0) {
        if (length(ini_val)>0){
          ini_val$B[,1] <<- ini_val$B %*% matrix(c(intercept, .colmeans_X), nc=1)
          ini_val$B[,-1] <<- ini_val$B[,-1] %*% diag(.colsds_X)
        }
        if (length(fit)>0){
          for (k in 1:length(fit)){
            fit[[k]]$B[,1] <<- fit[[k]]$B %*% matrix(c(intercept, .colmeans_X), nc=1)
            fit[[k]]$B[,-1] <<- fit[[k]]$B[,-1] %*% diag(.colsds_X)
          }
        }
        intercept <<- 1
        X <<- cbind(intercept,W)
      } else {
        X <<- W
        clear_all_results()
      }
      standardized_X <<- T
    },
    
    # fitting -----
    nllkh = function(L = NULL, R = NULL, B = NULL, Y_hat = NULL, with_r = NULL, limits = c(-Inf, Inf)){
      if (length(.obs)==0) .make_obs()
      
      if (is.null(Y_hat)) {
        if (is.null(L) || is.null(R) || (is.null(B) & p>0)){
          k <- which(r==with_r)
          L <- fit[[k]]$L
          R <- fit[[k]]$R
          if (p>0){
            B <- fit[[k]]$B
          }
        }
      } 
      val <- 0
      for (j in 1:dim[2]){
        i <- .obs$byj[[j]]
        if (length(i)==0) next
        if (is.null(Y_hat)) {
          y_hat <- L[i,,drop=F] %*% t(R[j,,drop=F])
          if (p>0) y_hat <- y_hat + X[i,fitting_cov,drop=F] %*% t(B[j,,drop=F])
        } else {
          y_hat <- Y_hat[i,j]
        }
        y <- Y[i,j]
        
        off_set <- 0
        if (scaled_Y) {
          multiplier <- if (length(.multiplier_Y)==1 ) .multiplier_Y else .multiplier_Y[j]
          off_set <- if (length(.off_set_Y)==1 ) .offset_Y else .offset_Y[j]
          y_hat <- y_hat*multiplier
          y <- y*multiplier
        }
        
        if (!is.infinite(limits[1])) y_hat[(y_hat + off_set)<limits[1]] <- limits[1] - off_set
        if (!is.infinite(limits[2])) y_hat[(y_hat + off_set)>limits[2]] <- limits[2] - off_set
        
        val <- val + sum((y - y_hat)^2)
      }
      return(val)
    },
    
    .make_initial = function(with_r=NULL, tol=.tol){
      if (is.null(with_r)) with_r <- max(r)
      if (length(ini_val)!=0 && ini_val$rank >= with_r){
        cat(sprintf('$ini_val already exists with r = %d... (%s)\n', ini_val$rank, Sys.time()))
        return(invisible(NULL))
      }
      .make_obs()
      estimate_pi()
      cat(sprintf('\nInitialization with rank %d ... (%s)\n', with_r, Sys.time()))
      
      if (p>0) {
        ans <- MC.covariate.iterative(Y = Y, X = X[,fitting_cov,drop=F], obs.index = .obs, pi.hat = pi_est$pi.hat,
                                      max.iter = 0, rank = with_r, iter.method = 'reg', miss.type = missing_model, tol = tol)
      } else {
        ans <- svd.missing(Z = Y, obs.index = .obs, rank = with_r, miss.type = missing_model, pi.hat = pi_est$pi.hat)
      }
      # ans$D <- NULL
      ans$Z <- NULL
      ans$converge <- NULL
      ini_val <<- ans
      cat(sprintf('$ini_val done! ... (%s)\n', Sys.time()))
    },
    
    .solve = function(with_rs, max_iter = .self$max_iter, tol = .tol){ 

      for (with_r in with_rs){
        k <- which(r==with_r)
        if (.solved[k]) next
        cat(sprintf('r = %02d, max_iter = %02d ... (%s)\n', with_r, max_iter,  Sys.time()))
        
        if (p!=0){

          fit[[k]] <<- MC.covariate.iterative(Y=Y, X=X[,fitting_cov,drop=F], obs.index=.obs, rank=with_r, pi.hat = pi_est$pi.hat, 
                                              B.ini=ini_val$B, L.ini=ini_val$L[,1:with_r, drop = F], R.ini=ini_val$R[,1:with_r, drop = F], 
                                              max.iter = max_iter, iter.method = 'reg', miss.type = missing_model, tol = tol
                                              )
        } else {
          fit[[k]] <<- lowrank.mc.iterative.reg(Z=Y, obs.index = .obs, 
                                                L.ini=ini_val$L[,1:with_r, drop = F], R.ini=ini_val$R[,1:with_r, drop = F], 
                                                max.iter = max_iter, tol = tol)
        }
        
        fit[[k]]$Z <<- NULL
        fit[[k]]$rmse <<- sqrt(nllkh(with_r = with_r)/obs_count)
        .solved[k] <<- T
      }
    },
    
    fitting = function(with_rs = r, max_iter = .self$max_iter, tol = .tol, force = F){
      if (length(.solved)==0) {
        .solved <<- rep(F, length(r)) 
        fit <<- list()
      }
      with_rs %<>% unique %>% sort
      
      new_r <- setdiff(with_rs, r)
      old_r <- setdiff(with_rs, new_r)
      
      fit_r <- new_r
      if (length(old_r)>0){
        if (force){
          index <- outer(r, old_r, '==') %>% apply(MARGIN = 2, which)
          clear_fit(index)
          new_r <- sort(c(new_r, old_r))
          fit_r <- new_r
        } else {
          fit_r <- c(fit_r, setdiff(old_r, r[.solved]))
        }
      }
      
      r <<- c(r, new_r)
      .solved <<- c(.solved, rep(F, length(new_r)))
      if (length(fit_r)==0){
        cat('\nNo new fitting needed\n')
        return(invisible(NULL))
      }
      .make_initial(with_r = max(fit_r))
      cat('\nfitting with rank(s):', fit_r, '\n' )
      .solve(with_rs = fit_r, max_iter = max_iter, tol = tol)
    },
    
    # post analysis (hypothesis testing, asyptotic dist., prediction) -----
    ht_beta = function(A_list, a_list = rep(list(0), length(A_list)), n_boots = 500, alpha = 0.05, with_r, track = T){
      if (p==0) return(invisible(NULL))
      k <- which(r==with_r)
      H <- MC.bs.ht(X = X[,fitting_cov,drop=F], Y = as.matrix(Y), obs.index = .obs, miss.type = missing_model,
                        B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                        A.list = A_list, a.list = a_list, n.boots = n_boots, alpha = alpha, track = track)
      H$rank <- with_r
      return(H)
    },
    
    ht_beta_simple = function(ps = NULL, n_boots = 500, alpha = 0.05, with_r, track = T){
      if (p==0) return(invisible(NULL))
      
      k <- which(r==with_r)
      H <- MC.bs.ht.simple(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                               B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                               alpha = alpha, n.boots = n_boots, track = track, selec = ps)
      if (is.null(ps)) ps <- 1:p
      H$describe <- sprintf('beta_j%d = 0, for all j', ps-(intercept!=0))
      H$rank <- with_r
      return(H)
    },
    
    ht_beta_allj = function(A.matrix = matrix(c(-1,1, rep(0, p - 2 + (intercept!=0))), nr=1), a.matrix = 0, n_boots = 500, alpha = 0.05, with_r, track = T){
      
      if (p==0) return(invisible(NULL))
      if (ncol(A.matrix)!=p) stop(sprintf('Incorrect column of matrix A , should be %d', p))
      
      pos <- A.matrix!=0
      describ <- character()
      for (l in 1:nrow(A.matrix)){
        u <- which(pos[l,])
        v <- A.matrix[l,u]
        if (length(v)==1) coeffs <- paste0(c('-', '')[(v>0) + 1], sprintf('%.2g', abs(v)))
        if (length(v)>1) coeffs <- paste0(c(c('-', '')[(v[1]>0) + 1], c(' - ', ' + ')[(v[-1]>0) + 1]), sprintf('%.2g', abs(v)))
        if(is.matrix(a.matrix)) {
          describ <- c(describ, paste0('H: ', paste0(coeffs, '*beta_j', u-(intercept!=0), collapse = ''), sprintf(' = a[%d,j]', l), ' for all j'))
        } else {
          if(length(a.matrix)==1) describ <- c(describ, paste0('H: ', paste0(coeffs, '*beta_j', u-(intercept!=0), collapse = ''), sprintf(' = %.4f', a.matrix), ' for all j'))
          if(length(a.matrix)==nrow(A.matrix)) describ <- c(describ, paste0('H: ', paste0(coeffs, '*beta_j', u-(intercept!=0), collapse = ''), sprintf(' = %.4f', a.matrix[l]), ' for all j'))
        }
      }
      
      k <- which(r==with_r)
      H <- MC.bs.ht.allj(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                             B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                             A.matrix = A.matrix, a.matrix = a.matrix, alpha = alpha, n.boots = n_boots, track = track)
      H$describe <- describ
      H$rank <- with_r
      return(H)
    },
    
    ht_beta_eachj = function(A_list, a_list = rep(list(0), length(A_list)), n_boots = 500, alpha = 0.05, with_r, track = T){
      if (p==0) return(invisible(NULL))
      k <- which(r==with_r)
      H <- MC.bs.ht.eachj(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                              B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                              A.list = A_list, a.list = a_list, n.boots = n_boots, alpha = alpha, track = track)
      H$describe <- 'A*beta_j = a'
      H$rank <- with_r
      return(H)
    },
    
    ht_beta_zero = function(n_boots = 500, alpha = 0.05, with_r, track = T){
      if (p==0) return(invisible(NULL))
      k <- which(r==with_r)
      H <- MC.bs.ht.all.zero(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                                 B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat, 
                                 n.boots = n_boots, alpha = alpha, track = track)
      H$describe <- 'beta_jp = 0 for all j,p'
      H$rank <- with_r
      return(H)
    },
    
    asymp_var = function(Zi = NULL, Zj = NULL, Bj = NULL, Bp = NULL, with_r){
      k <- which(r==with_r)
      var <- MC.asymp.cov.cond.LFX(X = X[,fitting_cov,drop=F], Y = Y, obs.index = .obs, miss.type = missing_model,
                                   B.est = fit[[k]]$B, L.est = fit[[k]]$L, R.est = fit[[k]]$R, pi.hat = pi_est$pi.hat,
                                   Zi = Zi, Zj = Zj, Bj = Bj, Bp = Bp)
      return(var)
    },
    
    predict = function(i = 1:.self$dim[1], j = 1:.self$dim[2], with_r = NULL, cross=T, limits = c(-Inf, Inf)){
      if (is.null(with_r)){
        print('please provide the rank. (use \'$r\' to find available ranks)')
        return(invisible(NULL))
      } else {
        k <- which(r==with_r)
      }
      
      if (cross){
        i <- unique(i) %>% sort
        j <- unique(j) %>% sort
      }
      
      if (p>0) {
        Xi <- X[i,fitting_cov,drop=F]
        Bj <- fit[[k]]$B[j,, drop = F]
      }
      
      
      Li <- fit[[k]]$L[i,, drop = F]
      Rj <- fit[[k]]$R[j,, drop = F]
      
      if (cross){
        ans <- Rj %*% t(Li)
        if (p>0) ans <- ans + Bj %*% t(Xi)
        if (scaled_Y)  ans <- ans*.multiplier_Y[j] + .offset_Y[j]
        ans <- t(ans)
        rownames(ans) <- i
        colnames(ans) <- j
      } else {
        ans <- if (with_r==1) Li*Rj else rowSums(Li*Rj)
        if (p>0) ans <- if(p==1) ans + Xi*Bj else ans + rowSums(Xi*Bj)
        if (scaled_Y){
          ans <- ans*.multiplier_Y[j] + .offset_Y[j]
        }
        names(ans) <- paste0( i, ',' , j)
      }
      if (!is.infinite(limits[1])) ans[ans<limits[1]] <- limits[1]
      if (!is.infinite(limits[2])) ans[ans>limits[2]] <- limits[2]
      return( ans)
    },
    
    # utilities -----
    summary = function(detail = F){
      cat('\nY')
      cat(sprintf('\n  %-16s: %d, %d', 'dim', .self$dim[1], .self$dim[2]))
      cat(sprintf('\n  %-16s: %s', 'scaled', ifelse(scaled_Y, 'Yes', 'No')))
      
      cat('\n\nX')
      if (ncol(X)==0) {
        cat(sprintf('\n  %-16s\n', 'No X matrix'))
        return(invisible(NULL))
      }
      cat(sprintf('\n  %-16s: %d, %d', 'dim', nrow(X), ncol(X)))
      cat(sprintf('\n  %-16s: %s', 'standardized', ifelse(standardized_X, 'Yes', 'No')))
      cat(sprintf('\n  %-16s: %s', 'intercept', ifelse(intercept!=0, round(intercept, 2), 'None')))
      cat(sprintf('\n  %-16s:', 'satats'))
      A <- apply(X, MARGIN = 2, quantile)
      rownames(A) <- c('Min', '1st', 'Med', '3rd', 'Max')
      cat(sprintf('\n  %-8s  %s', '', paste0(sprintf('%10s', colnames(A)), collapse = ' ')))
      for (i in 1:5) cat(sprintf('\n  %-8s  %s', rownames(A)[i], paste0(sprintf('%10.3f', A[i,]), collapse = ' ')))
      cat(sprintf('\n  %-8s  %s', 'Mean', paste(c(sprintf('%10s','')[intercept>0],sprintf('%10.3f', .colmeans_X)), collapse = ' ')))
      cat(sprintf('\n  %-8s  %s', 'SD', paste(c(sprintf('%10s','')[intercept>0],sprintf('%10.3f', .colsds_X)), collapse = ' ')))
      
      #cat(sprintf('\n  %-16s:\n', 'peek')); print(X[1:min(5, .self$dim[1]),,drop=F])
      
      cat('\n\nModel')
      cat(sprintf('\n  %-16s: %s', 'fitting with', paste0(colnames(X)[fitting_cov], collapse = ', ')))
      cat(sprintf('\n  %-16s: %s', 'missing model', missing_model))
      if (missing_model=='logistic') cat(sprintf('\n  %-16s: %s', 'pi fitting with', paste0(colnames(X)[pi_fitting_cov], collapse = ', ')))
      cat('\n')
    }
    # -----
  )
)

MC_alt_LS_rank_solver <- setRefClass('MC_alt_LS_rank_solver',
  contains = 'MC_alt_LS_solver',
  
  # ===== Methods =====
  methods = list(
    initialize = function(...){
      callSuper(...)
    },
    # rank estimation -----
    .get_r_bar = function(C = 1.5, ic_penalty_fun = NULL, tol=.tol, ...){
      estimate_pi()
      cat(sprintf('\ngetting r_bar for rank estimation ... (%s)\n', Sys.time()))
      if (p>0) {
        ans <- MC.covariate.iterative(Y = Y, X = X[,fitting_cov,drop=F], obs.index = .obs, pi.hat = pi_est$pi.hat,
                                      max.iter = 0, rank = 0, iter.method = 'reg', miss.type = missing_model, tol = tol)
      } else {
        ans <- svd.missing(Z = Y, obs.index = .obs, rank = 0, miss.type = missing_model, pi.hat = pi_est$pi.hat)
      }
      
      if (missing_model == 'MCR'){
        alp <- .self$pi_est$pi.hat
      }
      
      if (missing_model == 'logistic'){
        alp <- exp(.self$pi_est$gam.hat['(Intercept)'])
      }
      
      ic <- IC(n = .self$dim[1], m = .self$dim[2], d = ans$D, pfun = ic_penalty_fun, alpha_n=alp, ...)
      d_ic <- diff(ic)
      rank_est$r_bar <<- max(as.integer(C*(match(T, d_ic>0))), 5)
      cat(sprintf('$rank_est$r_bar = %d ... (%s)\n', rank_est$r_bar, Sys.time()))
    },
    
    .solve_rank = function(r_bar, max_iter = .self$max_iter, tol = .tol){ 
      
      if (p>0){
        W <- Y - X[,fitting_cov,drop=F] %*% t(ini_val$B)
        W[-.obs$all] <- 0
        W <- Matrix(W, sparse = TRUE)
      }
      for (with_r in 1:r_bar){
        if (!is.na(rank_est$mse[with_r])) next
        cat(sprintf('r = %02d, max_iter = %02d ... (%s)\n', with_r, max_iter,  Sys.time()))
        if (p>0){
          ans <- lowrank.mc.iterative.reg(Z=W, obs.index = .obs, 
                                          L.ini=ini_val$L[,1:with_r, drop = F], 
                                          max.iter = max_iter, tol = tol)
        } else {
          ans <- lowrank.mc.iterative.reg(Z=Y, obs.index = .obs, 
                                          L.ini=ini_val$L[,1:with_r, drop = F], 
                                          max.iter = max_iter, tol = tol)
          
        }
        rank_est$mse[with_r] <<- nllkh(L = ans$L, R = ans$R, B = ini_val$B)/obs_count
      }
    },
    
    rank_estimation = function(r_bar = NULL, C_rbar = 1.5, penalty_fun = NULL, penalty_names = NULL, ...){
      
      if (is.null(r_bar) || r_bar<1){
        if (is.null(rank_est$r_bar)) .get_r_bar(C = C_rbar)
      } else {
        rank_est$r_bar <<- as.integer(r_bar)
      }
      
      len_pen_old <- 0
      r_bar_old <- 0
      if (!is.null(rank_est$eIC)){
        len_pen_old <- ncol(rank_est$eIC)
        r_bar_old <- nrow(rank_est$eIC)
        rank_est$r_bar <<- max(rank_est$r_bar, r_bar_old)
        if (rank_est$r_bar>r_bar_old) {
          rank_est$eIC <<- rbind(rank_est$eIC, matrix(nr=rank_est$r_bar-r_bar_old, nc=len_pen_old))
        }
      }
      
      rank_est$mse <<- c(rank_est$mse, rep(NA, rank_est$r_bar - length(rank_est$mse)))
      
      len_pen <- length(eIC(n = 1, m = 1, mse = 1, r = 1, pfun = penalty_fun, ...))
      n_old <- 0
      for (name in penalty_names){
        if (name %in% colnames(rank_est$eIC)) n_old <- n_old+1
      }
      
      if (length(penalty_names)<len_pen){
        ext <- 1:(len_pen - length(penalty_names)) + length(penalty_names)- n_old + len_pen_old
        penalty_names <- c(penalty_names, paste0('pen[', ext, ']'))
      }
      penalty_names <- penalty_names[1:len_pen]
      
      .make_initial(with_r = rank_est$r_bar) 
      cat(sprintf('\nRank estimation with max rank = %d ... (%s)\n', rank_est$r_bar, Sys.time()))
      .solve_rank(r_bar=rank_est$r_bar)
      
      if (missing_model == 'MCR'){
        alp <- .self$pi_est$pi.hat
      }
      
      if (missing_model == 'logistic'){
        alp <- exp(.self$pi_est$gam.hat['(Intercept)'])
      }
      
      cat(sprintf('\nCalculating eIC with penalty %s ... (%s)\n', paste(penalty_names, collapse = ', '), Sys.time()))
      eic <- eIC(n = .self$dim[1], m = .self$dim[2], mse = rank_est$mse, pfun = penalty_fun, alpha_n=alp, ...)
      colnames(eic) <- penalty_names
      
      for (name in penalty_names){
        if (name %in% colnames(rank_est$eIC)){
          rank_est$eIC[,name] <<- eic[,name, drop=F]
        } else {
          rank_est$eIC <<- cbind(rank_est$eIC,eic[,name, drop=F])
        }
      }
      row.names(rank_est$eIC) <<- paste0('r=', 1:rank_est$r_bar)
      rank_est$est <<- apply(rank_est$eIC, MARGIN = 2, which.min)
      cat(sprintf('$rank_est with max rank = %d done! ... (%s)\n', rank_est$r_bar, Sys.time()))
    }
    #-----
  )
)
  