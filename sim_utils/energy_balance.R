energy_balance<-function(y,treat,X,taus,stable=FALSE,lambda=NULL,normalize=FALSE,min.w=1e-8,
                         improved=FALSE,moments=FALSE,geq0=TRUE,sum1=TRUE){
  
  n<- length(treat)
  s.weights<-rep(1,n)
  y_grid<-sort(unique(y))
  
  if(normalize==TRUE) X<-scale(X,center = TRUE,scale = TRUE)
  Dmat <- as.matrix(dist(X))        # Euclidean distance matrix
  
  # One n*2 logical matrix for treat
  levels_treat <- as.numeric(levels(factor(treat))) 
  tmat <- vapply(levels_treat, function(t) treat == t, logical(n))
  nt <- colSums(tmat)     # num of treat=0 and treat=1
  
  
  # Two vector: (1/n0,....,1/n0,0,...,0);(0,...,0,1/n1,....,1/n1)
  w_diffsample<- setNames(lapply(levels_treat, function(t) s.weights*tmat[,t+1]/nt[t+1]), levels_treat) 
  # One vector: (1/n,...,1/n)
  w_allsample<- as.matrix(s.weights/n)    
  w_diffsum <- w_diffsample[[1]]+w_diffsample[[2]]
  
  # Two matrix: 
  M2_array <- vapply(levels_treat, function(t)  -2*tcrossprod(w_diffsample[[(t+1)]]) * Dmat, diag(n))
  #M2_array2 <- w_diffsample[[1]]%*%t(w_diffsample[[2]])* Dmat
  M1_array <- vapply(levels_treat, function(t) 2 * w_diffsample[[t+1]] * Dmat %*% w_allsample, w_allsample)
  
  M2 <- rowSums(M2_array, dims = 2)
  M1 <- rowSums(M1_array)
  
  
  if (improved==TRUE) {
    #all_pairs <- combn(levels_treat, 2, simplify = FALSE)
    M2_pairs_array <-  -2 * tcrossprod(w_diffsample[[1]]-w_diffsample[[2]]) * Dmat
    M2 <- M2 + M2_pairs_array 
  }
  
  
  if (stable==TRUE) diag(M2) <- diag(M2) + lambda/ n
  
  Amat <- rbind(diag(n))
  #geq0='FALSE';sum1='FALSE'
  if(geq0==TRUE){
    min.w = 1e-8
    lvec <- rep(min.w, n)
    uvec <- rep(Inf,n)
  }else{
    lvec <- rep(-Inf, n)
    uvec <- rep(Inf,n)
  }
  
  if(sum1==TRUE){
    Amat<-rbind(Amat,t(s.weights * tmat))
    lvec<-c(lvec,nt)
    uvec<-c(uvec,nt)
  }
  
  if(moments==TRUE){
    Amat <- rbind(Amat,t(X))
    lvec <- c(lvec,colMeans(X))
    uvec <- c(uvec,colMeans(X))
  }
  
  settings <- osqpSettings(verbose=FALSE,max_iter = 500)
  opt.out <- solve_osqp(P = M2, q = M1, A = Amat, l = lvec, u = uvec,settings)
  
  w <- opt.out$x
  converge_status<-opt.out$info$status
  
  ATE <- get_ate(y,treat,w)
  cdf_qte <- get_qte(y,treat,w,taus,y_grid)
  QTE <- cdf_qte$QTE
  Q0 <- cdf_qte$Quantiles0
  Q1 <- cdf_qte$Quantiles1
  
  F0 <- cdf_qte$F0
  F1 <- cdf_qte$F1
  
  results<-list(ATE=ATE,QTE=QTE,Q0=Q0,Q1=Q1,F0 = F0,F1 = F1,w=w,taus=taus,converge_status=converge_status)
  return(results)
}