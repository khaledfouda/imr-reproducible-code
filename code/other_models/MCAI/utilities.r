library(MASS)

# rank estimation for MC problem
#-------------------------------------------------------------
ER <- function(d = c(1,1), r.min = 1, r.max = length(d)-1){
  return(d[r.min:r.max]/d[r.min:r.max+1])
}

IC <- function(n = 2, m = 2, d = c(1,1), r.min = 1, r.max = length(d), pfun = NULL, ...){ 
  d.sq <- d^2
  if (is.null(pfun)){
    pfun <- function(n, m, alpha_n=1, delta_h = 0.2, C_h = 1, ...){
      u <- sqrt(alpha_n)*m*n/(m+n)
      pen <- c(outer(C_h, n^delta_h, '*'))/u
      return(pen)
    }
  }
  pen <- pfun(n,m, ...)
  ic <- log( sum(d.sq) - cumsum(d.sq[1:r.max])[r.min:r.max] ) - log(m*n) +  outer((r.min:r.max), pen, '*')
  rownames(ic)<- r.min:r.max
  return(ic)
}

eIC <- function(n = 2, m = 2, mse = 1, r = 1:length(mse), pfun = NULL, ...){
  if (is.null(pfun)){
    pfun <- function(n, m, alpha_n=1, obr=1, delta_h = 0.1, C_h = 0.9, ...){
      u <- alpha_n*m*n/(m+n)
      pen <- c(outer(C_h, n^delta_h, '*'))/sqrt(u)
      return(pen)
    }
  }
  pen <- pfun(n,m,...)

  eic <- log(mse) + outer(r,pen, '*')
  rownames(eic)<- r
  return(eic)
}
#-------------------------------------------------------------

# solving missing value problem with linear model
#-------------------------------------------------------------
lowrank.mc.iterative.svd <- function(Z, rank, miss = NULL, obs.index = NULL, pi.hat = NULL, miss.type = c('MCR', 'logistic')[1], X = NULL,
                                     Z.ini = NULL, tol = 1e-6, max.iter = 50){
  
  # low rank approximation via iterative SVD
  find.missing <- function(){
    if (is.null(miss)) {
      reurn(Z.est[-obs.index$all])
    } else {
      return(Z.est[miss])
    }
  }
  
  n <- dim(Z)[1]
  m <- dim(Z)[2]
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  
  Z.est <- Z.ini
  err <- tol + 1
  iter <- 0
  ans <- list(Z = Z.est, L = NULL, R = NULL)
  while (err>tol & iter<max.iter){
    if(is.null(Z.est)) {
      Z.comp <- NULL
    } else {
      Z.comp <- find.missing()
    }
    ans <- svd.missing(Z = Z, rank = rank, miss = miss, obs.index = obs.index, miss.type = miss.type, pi.hat = pi.hat, Z.comp = Z.comp)
    Z.new <- ans$L %*% t(ans$R)
    err <- ifelse(is.null(Z.est), tol + 1, mean((Z.est - Z.new)^2))
    Z.est <- Z.new
    iter <- iter + 1
  }
  
  return(list(Z = Z.est, L = ans$L, R = ans$R, rank = rank, 
              converge = list(iter = iter, reach.max = iter==max.iter, err = err)
  ))
}

lowrank.mc.iterative.reg <- function(Z, miss = NULL, obs.index = NULL, L.ini, R.ini = NULL, 
                                     tol = 1e-6, max.iter = 50, epsl=1e-16){
  
  # low rank approximation via iterative regression
  
  L <- L.ini
  R <- R.ini
  Z.est <- if(is.null(R)) NULL  else L %*% t(R)
  iter <- 0
  err <- tol + 1
  while (iter<max.iter && err>tol){
    R <- reg.missing(Y = Z, X = L, miss = miss, obs.index_byj = obs.index$byj, epsl = epsl)
    L <- reg.missing(Y = t(Z), X = R, miss = if(is.null(miss)) NULL else t(miss), obs.index_byj = obs.index$byi, epsl = epsl)
    
    Z.new <- L %*% t(R)
    if (is.null(Z.est)){
      err <- tol + 1
    } else {
      err <- if (!is.null(miss)) mean((Z.est[!miss] - Z.new[!miss])^2) else mean((Z.est[obs.index$all] - Z.new[obs.index$all])^2)
    }
    Z.est <- Z.new
    iter <- iter + 1
  }
  return(list(Z = Z.est, L = L, R = R, rank = dim(L)[2], 
              converge = list(iter = iter, reach.max = iter==max.iter, err = err)
  ))
}

pi.est <- function(miss = NULL, obs.index = NULL, miss.type = c('MCR', 'logistic')[1], X = NULL, detail=F){
  
  if (!is.null(X)){
    vars <- apply(X, MARGIN = 2, var)
    X <- X[,vars!=0]
  }
  
  n <- if (!is.null(miss)) nrow(miss) else length(obs.index$byi)
  m <- if (!is.null(miss)) ncol(miss) else length(obs.index$byj)
  
  if (miss.type == 'MCR'){
    
    pi.hat <- if (!is.null(miss)) mean(!miss) else sum(lengths(obs.index$byi))/(n*m)
    return(list(pi.hat = pi.hat))
    
  } else if (miss.type == 'logistic'){
    if (is.null(X) || ncol(X)==0) stop('No X for logistic model')
    n.obs <- if (!is.null(miss)) rowSums(!miss) else lengths(obs.index$byi)
    n.obs <- cbind(n.obs,m-n.obs)
    
    fit <- glm(n.obs~1+X, binomial(link = "logit"))
    res <- list(pi.hat = fit$fitted.values, gam.hat = fit$coefficients)
    if (detail) res$fit <- fit
    return(res)
  } else {
    stop('Invalid missing_type')
  }
}

svd.missing <- function(Z, miss = NULL, obs.index = NULL, Z.comp = NULL, rank, 
                        pi.hat = NULL, miss.type = c('MCR', 'logistic')[1], X = NULL){
  
  n <- dim(Z)[1]
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  
  
  if (miss.type=='MCR' && pi.hat==1) {
    svd.result <- svd(Z, nu = rank, nv= rank)
    return(
      list( L = svd.result$u * sqrt(n), 
            R = svd.result$v %*% diag(x = svd.result$d[1:rank], nrow = rank) / sqrt(n), 
            D = svd.result$d,
            rank = rank)
    )
  }
  
  # we pad the provided value, or we pad 0 and divide the matrix by pi-hat
  if (!is.null(Z.comp)){
    if (is.null(miss)) {
      Z[-obs.index$all] <- Z.comp
    } else { Z[miss] <-  Z.comp }
  } else {
    if (is.null(miss)){
      Z[-obs.index$all] <- 0
      Z <- Z/pi.hat
    } else {
      Z[miss] <- 0
      Z <- Z/pi.hat
    }
  }
  
  svd.result <- svd(Z, nu = rank, nv= rank)
  if (rank==0){
    return(
      list( D = svd.result$d,
            rank = rank)
    )
  }
  return(
    list( L = svd.result$u * sqrt(n), 
          R = svd.result$v %*% diag(x = svd.result$d[1:rank], nrow = rank) / sqrt(n), 
          D = svd.result$d,
          rank = rank)
  )
}

reg.missing <- function(Y, X, miss = NULL, obs.index_byj = NULL, penalty = 0, epsl = 1e-16){
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]

  q <- if (!is.null(miss)) mean(!miss) else sum(lengths(obs.index_byj))/(n*m)
  
  if (q ==1) return( t(ginv(t(X) %*% X) %*% (t(X) %*% Y)) )
  
  p <- dim(X)[2]
  
  # solve each beta_j by using only those subject with non-missing response on j-th item

  B <- matrix(0, nr=m,nc=p)
  for (j in 1:m){
    if (is.null(miss)) {
      i<- obs.index_byj[[j]]
      k <- length(i)
    } else {
      i <- !miss[,j]
      k <- sum(i)
    }
    if (k==0) next
    Xr <- X[i, , drop = F]
    Yr <- Y[i,j, drop = F]
    H <- t(Xr) %*% Xr/k
    if (penalty>0) H <- H + diag(penalty, p)
    b <- as.vector( ginv(H) %*% (t(Xr) %*% Yr)/k )
    b[abs(b)<epsl] <- 0
    B[j,] <- b
  }
  
  return(B)
}
