# The argument "miss" is the logic matrix that indicates the missing position
# obs.index is a list contains tow long lists obs.index$byi and obs.index$byj st obs.index$byi[[i]] or obs.index$byj[[j]] is a vector of observed index

# methods for MC problem
#-------------------------------------------------------------
MC.complete <- function(Y, X = NULL, rank = NULL, return.complete.hidden = T){
  ans <- list()
  n <- nrow(Y)
  if (!is.null(X)) {
    B <- ginv(t(X) %*% X, tol = 1e-6) %*% (t(X) %*% Y)
    Y <- Y - X %*% B
    ans$B <- t(B)
  }
  svd.result <- svd(Y, nu = rank, nv= rank)
  L = svd.result$u * sqrt(n)
  R = svd.result$v %*% diag(x = svd.result$d[1:rank], nrow = rank) / sqrt(n)
  if (return.complete.hidden) ans$Z <- L %*% t(R)
  return( c(ans, 
            list(L = L, R = R, D = svd.result$d, rank = rank, converge =list( iter = 0, reach.max = NA, err = NA) ) 
            )
  )
}



MC.covariate.iterative <- function(Y, X, miss = NULL, obs.index = NULL, rank = NULL,
                                   B.ini = NULL, L.ini = NULL, R.ini = NULL,
                                   max.iter = 50, tol = 1e-6, iter.method = c('reg', 'svd')[1],
                                   pi.hat = NULL, miss.type = c('MCR', 'logistic')[1]){
  
  # alternatively solving B and Z until converge
  # solving Z via regression or SVD
  
  # the given missing information could be a full matrix, 'miss', with (T,F) missing indicator 
  # or a list of observed index, obs.index, including indices based on row ($byi), based on cloumn ($byj) and based on the matrix ($all) 
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  
  if (miss.type=='MCR' && pi.hat==1) return(MC.complete(Y = Y, X = X, rank = rank))
  
  residul.matrix <- function(X=NULL,B=NULL,Z=NULL){
    W <- if (is.null(Z)) Y-X %*% t(B) else Y-Z
    if (!is.null(miss)){
      W[miss] <- 0
    } else {
      W[-obs.index$all] <- 0
    }
    if (inherits(Y, 'sparseMatrix')) W <- Matrix(W, sparse = TRUE)
    return(W)
  }
  
  if (is.null(B.ini) && (is.null(L.ini) || is.null(R.ini))) B.ini <- reg.missing(Y, X, miss,  obs.index_byj = obs.index$byj)
  
  ans <- list(L=L.ini, R=R.ini)
  B <- B.ini
  
  if (is.null(L.ini) || is.null(R.ini)){
    W <- residul.matrix(X = X , B = B.ini)
    if (iter.method=='reg' && (!is.null(L.ini) && !is.null(R.ini))){
      if (!is.null(L.ini)) ans$R <- reg.missing(W, L.ini, miss = miss, obs.index_byj = obs.index$byj)
      if (!is.null(R.ini)) ans$L <- reg.missing(t(W), R.ini, miss = if (is.null(miss)) NULL else t(miss), obs.index_byj = obs.index$byi)
    }
    if ( (is.null(L.ini) && is.null(R.ini)) || iter.method=='svd') {
      ans <- svd.missing(Z = W, miss = miss, obs.index = obs.index, rank = rank, pi.hat = pi.hat, miss.type = miss.type)
    }
  }
  
  if (is.null(B.ini)) B <- 0
  
  if (rank==0){
    ans$B <- B
    return(ans)
  }
  ans$Z <- ans$L %*% t(ans$R)
  Yh <- X%*%t(B) + ans$Z
  err.B <- NA
  err.Z <- NA
  err <- tol + 1
  iter <- 0
  while (err>tol && iter<max.iter) {
    W <- residul.matrix(Z = ans$Z)
    B <- reg.missing(Y = W, X = X, miss = miss, obs.index_byj = obs.index$byj)
    W <- residul.matrix(X = X, B = B)
    if (iter.method == 'reg') ans <- lowrank.mc.iterative.reg(Z = W, miss = miss, obs.index = obs.index, L.ini = ans$L, R.ini = ans$R, tol = tol, max.iter = 1)
    if (iter.method == 'svd') ans <- lowrank.mc.iterative.svd(Z = W, miss = miss, obs.index = obs.index, miss.type = miss.type, Z.ini = ans$Z, tol = tol, max.iter = 1, rank = rank)
    Yh.new <- X%*%t(B) + ans$Z
    err <- if(!is.null(miss)) mean((Yh[!miss] - Yh.new[!miss])^2) else mean((Yh[obs.index$all] - Yh.new[obs.index$all])^2)
    Yh <- Yh.new
    iter <- iter + 1
  }
  ans$B <- B
  ans$converge <- list( iter = iter, reach.max = (iter==max.iter), err = err)
  return(ans)
}

#-------------------------------------------------------------

# Hypothesis testing
#-------------------------------------------------------------
MC.bs.ht <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                         A.list, a.list, n.boots = 1000, alpha = 0.05, track = T){
  
  # Hypothesis testing on 
  #
  #  H0: A*vec(t(B))=a 
  #
  # w/ normal multiplier bootstrap method 
  # Do the above test for given pairs A=A.list[i], a=a.list[i] i=1, 2, 3,...
  
  pivot <- function(A,a,B){ max(abs(A %*% matrix(t(B), ncol = 1) - a)) }
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  
  r <- dim(L.est)[2]
  p <- dim(B.est)[2]
  
  n.test <- length(A.list)
  
  test.names <- paste0('test_', 1:n.test)
  test.stat <- numeric(n.test)
  for (l in 1:n.test){
    test.stat[l] <- pivot(A.list[[l]], a.list[[l]], B.est)
  }
  bs.stat <- array(dim=c(n.test, n.boots), dimnames = list(test.names, 1:n.boots))
  pv <- numeric(n.test)
  names(pv) <- test.names
  names(test.stat) <- test.names
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))
  
  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)
  
  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est
  
  for (k in 1:n.boots){
    for (l in 1:n.test){
      tmp <- Sig.X.inv.Xt %*% (rnorm(n)*W)/n
      bs.stat[l,k] <- pivot(A.list[[l]], 0, t(tmp))
    }
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- rowMeans(bs.stat > test.stat)
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat, A = A.list, a = a.list)
  )# 'FALSE' means reject null
}

MC.bs.ht.eachj <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                               A.list, a.list, n.boots = 1000, alpha = 0.05, track = T){
  
  # A: q-by-p matrix, a: q-by-m vector
  # Hypothesis testing on 
  #
  #   H0_j: A*b_j=a 
  #
  # for each b_j=t(B[j,]) w/ normal multiplier bootstrap method 
  # Do the above test for given pairs A=A.list[i], a=a.list[i] i=1, 2, 3,...
  
  pivot <- function(A,a,B){ apply(abs(A %*% t(B) - c(a)), MARGIN = 2, max)  }
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  
  r <- dim(L.est)[2]
  p <- dim(B.est)[2]
  
  n.test <- length(A.list)
  test.names <- list(paste0('beta_',1:m), paste0('test_',1:n.test))
  test.stat <- array(dim=c(m, n.test), dimnames = test.names)
  for (l in 1:n.test){
    test.stat[,l] <- pivot(A.list[[l]], a.list[[l]], B.est)
  }
  pv <- array(dim=c(m, n.test), dimnames = test.names)
  bs.stat <- array(dim=c(m, n.test, n.boots), dimnames = c(test.names, list(1:n.boots)))
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))
  
  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)
  
  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est
  
  for (k in 1:n.boots){
    for (l in 1:n.test){
      tmp <- Sig.X.inv.Xt %*% (rnorm(n*m)*W)/n
      bs.stat[,l,k] <- pivot(A.list[[l]], 0, t(tmp))
    }
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- apply(bs.stat>c(test.stat), MARGIN = c(1,2), FUN = mean)
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat, A = A.list, a = a.list)
  )# 'FALSE' means reject null
}

MC.bs.ht.allj <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                              A.matrix, a.matrix = 0, n.boots = 1000, alpha = 0.05, track = T){
  
  # A.matrix: n.test x p
  # a.matrix: n.test x m 
  #
  # H0: <v, b_j>=a_j for all j
  #
  # where <,> is the inner-product of vectors, v=A.matrix[k,], a_j=a[k,j] for k = 1, 2, ..., n.test.
  #
  # Best use for testing contrasts, combining groups.
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  
  r <- dim(L.est)[2]
  p <- dim(B.est)[2]
  
  n.test <- nrow(A.matrix)
  
  test.names <- paste0('test_', 1:n.test)
  test.stat <- apply( abs(A.matrix %*% t(B.est) - a.matrix), MARGIN = 1, max)
  bs.stat <- array(dim=c(n.test, n.boots), dimnames = list(test.names, 1:n.boots))
  pv <- numeric(n.test)
  names(pv) <- test.names
  names(test.stat) <- test.names
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  A.Sig.X.inv.Xt <- A.matrix %*% solve(t(pi.hat*X)%*%X/n, t(X))
  
  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)
  
  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est
  
  N <- prod(dim(A.Sig.X.inv.Xt))
  for (k in 1:n.boots){
    tmp <- abs((rnorm(N)*A.Sig.X.inv.Xt) %*% W/n)
    bs.stat[,k] <- apply(tmp, MARGIN = 1, max)
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- rowMeans(bs.stat>test.stat)
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat, A = A.matrix, a = a.matrix)
  )# 'FALSE' means reject null
}

MC.bs.ht.simple <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                                n.boots = 1000, alpha = 0.05, track = T, selec=NULL){
  
  # For each p in the 'selec', test H0: B[, p]=0 
  # Use MC.bs.ht() for more general tests
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  
  r <- dim(L.est)[2]
  p <- dim(B.est)[2]
  
  if(is.null(selec)) selec <- (1:p)
  selec <- sort(unique(selec[selec<=p]))
  p <- length(selec) # this is now the number of vectors of B in the tests
  
  
  test.names <- sprintf('beta_%d', selec)
  
  test.stat <- apply(abs(B.est[, selec, drop=F]), MARGIN = 2, max)
  bs.stat <- array(dim=c(p, n.boots), dimnames = list(test.names, 1:n.boots))
  pv <- numeric(p); 
  names(pv) <- test.names
  names(test.stat) <- test.names
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))
  Sig.X.inv.Xt <- Sig.X.inv.Xt[selec,]
  
  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)
  
  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est
  
  for (k in 1:n.boots){
    tmp <- (Sig.X.inv.Xt*rnorm(p*n)) %*% W/n
    bs.stat[,k] <- apply(abs(tmp), MARGIN = 1, max)
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- rowMeans(bs.stat>test.stat)
  
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat)
  )# 'FALSE' means reject null
}

MC.bs.ht.all.zero <- function(X, Y, miss = NULL, obs.index=NULL, miss.type = c('MCR', 'logistic')[1], B.est, L.est, R.est, pi.hat=NULL,
                                  n.boots = 1000, alpha = 0.05, track = T){
  
  # A simple testing on H0: B=0. This means all the covariates have no effect
  # Use MC.bs.ht() for more general tests
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  
  r <- dim(L.est)[2]
  p <- dim(B.est)[2]
  
  test.stat <- max(abs(B.est))
  bs.stat <- numeric(n.boots)
  pv <- numeric(1)
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  Sig.X.inv.Xt <- solve(t(pi.hat*X)%*%X/n, t(X))
  
  Z.est <- L.est %*% t(R.est)
  W <- (Y - X %*% t(B.est) - Z.est)
  
  if (is.null(miss)){
    W[-obs.index$all] <- 0
  } else {
    W[miss]  <- 0
  }
  W <- W + pi.hat*Z.est
  
  pv <- 0
  for (k in 1:n.boots){
    tmp <- Sig.X.inv.Xt %*% (rnorm(n)*W/n)
    bs.stat[k] <- max(abs(tmp))
    if (track){
      if (mod(k,10)==0) cat('-')
      if (mod(k,200)==0) cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
    }
  }
  if (track & (mod(k,200)!=0)) {
    cat(sprintf('%s', rep('x', 20 - n.boots%%200%/%10)), sep = '')
    cat(sprintf(' %d/%d bootstraps done (%s)\n', k, n.boots, Sys.time()))
  }
  pv <- mean(bs.stat>test.stat)
  
  return(
    list(h.null = pv>alpha, pv = pv, test.stat = test.stat, bs.stat = bs.stat)
  )# 'FALSE' means reject null
}

#-------------------------------------------------------------

# Asymptotic covariance matrix
#-------------------------------------------------------------
MC.asymp.cov.cond.LFX <- function(X, Y, miss = NULL, obs.index = NULL, miss.type = c('MCR', 'logistic')[1],
                                  B.est, L.est, R.est, pi.hat=NULL, Zi = NULL, Zj = NULL, Bj = NULL, Bp = NULL){
  
  # input:
  # The data, the estimators and the desired index vectors: Zi, Zj, Bj, Bp. e.g. Zi=c(1,2,3), Zj=c(1,3,5), Bj=c(1:10), Bp=c(1:3)
  # Usually Bp is 1:p (or leave it NULL), so you'll get the full cov matrix of B_j
  
  # output:
  # the estimated variance of the estimator Z_ij for any i in Zi and j in Zj will be in out$Z.var
  # the estimated variance of the estimator EY_ij for any i in Zi and j in Zj will be in out$EY.var
  # the estimated covariance matrix of the estimator B_j[Bp] for any j in Bj will be in out$B.cov
  
  p <- dim(X)[2]
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  r <- dim(L.est)[2]
  
  # make sure each index is valid and unique
  Zi <- sort(unique(Zi[Zi<=n]))
  Zj <- sort(unique(Zj[Zj<=m]))
  Bj <- sort(unique(Bj[Bj<=m]))
  Bp <- sort(unique(Bp[Bp<=p]))
  
  # count how many we need, and get all if null
  zn <- ifelse(is.null(Zi), n, length(Zi))
  zm <- ifelse(is.null(Zj), m, length(Zj))
  bm <- ifelse(is.null(Bj), m, length(Bj))
  bp <- ifelse(is.null(Bp), p, length(Bp))
  
  if (bm*bp==0){
    B.cov <- NULL
    bm <- 0
  } else {
    B.cov <- array(0, dim = c(bp, bp, bm))
    p.names <- if (bp<p) Bp else 1:bp
    j.names <- if (bm<m) Bj else 1:bm
    dimnames(B.cov) <- list(p.names, p.names, j.names)
  }
  
  if (zn*zm==0){
    Z.var <- NULL
    EY.var <- NULL
    zm <- 0
  } else {
    Z.var <- matrix(0, nr=zn, nc=zm)
    i.names <- if (zn<n) Zi else 1:zn
    j.names <- if (zm<m) Zj else 1:zm
    dimnames(Z.var) <- list(i.names, j.names)
  }
  
  if (is.null(pi.hat)) pi.hat <- pi.est(miss = miss, obs.index = obs.index, miss.type = miss.type, X = X)$pi.hat
  
  Z.est <- L.est%*%t(R.est)
  
  # error variance estimation
  if (!is.null(miss)){
    W <- Y - X %*% t(B.est)
    W[miss] <- Z.est[miss]
  } else {
    W <- Z.est
    U <- X %*% t(B.est)
    W[obs.index$all] <- Y[obs.index$all] - U[obs.index$all]
  }
  
  d <- svd(W, nu = 0, nv = 0)$d^2
  df <- m*n*mean(pi.hat) - m*p - (m+n-r)*r - 2*sum( d[(r+1):min(n,m)]/outer(-d[(r+1):min(n,m)], d[1:r], '+')  )
  sigE <- sum(d[(r+1):min(n,m)])/df
  rm(d, W, df)
  
  # constructing var and cov matrix
  sigX.inv <- ginv(t(X*pi.hat)%*%X/n)
  sigR.inv <- ginv(t(R.est) %*% R.est/m)
  sigL.inv <- ginv(t(L.est*pi.hat) %*% L.est/n)
  
  js <- if (is.null(Bj) || is.null(Zj)) 1:m else sort(unique(c(Bj, Zj)))
  
  zj <- 1
  bj <- 1
  for (j in js){
    tmp <-  sigX.inv %*% ( t(X) %*% (X * (pi.hat*Z.est[,j])^2 )) %*% sigX.inv/n #middle part of zeta_ij
    if (bj <= bm && (is.null(Bj) || Bj[bj]==j)) {
      B.cov[,,bj] <- if (is.null(Bp)) (sigE*sigX.inv + tmp)/n else (sigE*sigX.inv[Bp, Bp] + tmp[Bp, Bp])/n
      bj <- bj +1
    }
    if (zj <= zm && (is.null(Zj) || Zj[zj]==j)){
      Z.var[,zj] <- if (is.null(Zi)) rowSums( (X %*% tmp) * X )/n else rowSums( (X[Zi,,drop=F] %*% tmp) * X[Zi,,drop=F] )/n
      zj <- zj +1
    } 
  }
  if (zn*zm==0)  return(list( B.cov = B.cov , Z.var = Z.var, EY.var = EY.var))
  if (!is.null(Zi)) {
    L.est <- L.est[Zi,,drop=F]
    X <- X[Zi,,drop=F]
    if (miss.type == 'MCR') {
      pi.hat <- rep(pi.hat, length(Zi))
    } else {
      pi.hat <- pi.hat[Zi]
    }
  }
  if (!is.null(Zj)) R.est <- R.est[Zj,,drop=F]
  tmp <- sigE*(outer( 1/pi.hat, rowSums((R.est %*% sigR.inv) * R.est )/m, '*') + rowSums((L.est %*% sigL.inv) * L.est )/n)
  Z.var <- Z.var + tmp
  EY.var <- tmp + sigE * rowSums( (X %*% sigX.inv) * X )/n
  dimnames(EY.var) <- dimnames(Z.var)
  return(list( B.cov = B.cov , Z.var = Z.var, EY.var = EY.var))
}

#-------------------------------------------------------------