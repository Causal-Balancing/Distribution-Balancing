library(WeightIt)

#Treat=t
SLE<-function(y,Treat,X,d,taus,y_grid){
  
  X<-as.matrix(X)
  
  if(missing(y_grid)){y_grid<-sort(unique(y))}
  RX<-powers_of_x(X,d)
  pscore.reg <- glm(Treat ~ RX-1,family=binomial)
  pscore <- fitted(pscore.reg)
  

  weights0 <- (1-D) / (1-pscore)
  weights0 <- weights0/sum(weights0)
  weights1 <- D / (pscore) 
  weights1 <- weights1/sum(weights1)
  weights <- weights0+weights1
  
  
  
  ATE <- get_ate(y,Treat,weights)
  cdf_qte <- get_qte(y,Treat,weights,taus,y_grid)
  QTE <- cdf_qte$QTE
  Q0 <- cdf_qte$Quantiles0
  Q1 <- cdf_qte$Quantiles1
  
  F0 <- cdf_qte$F0
  F1 <- cdf_qte$F1
  
  results<-list(ATE=ATE,QTE=QTE,Q0=Q0,Q1=Q1,F0 = F0,F1 = F1,weights=weights,taus=taus)
  return(results)
}
  


powers_of_x<-function(X,d){
  
  X<-as.matrix(X)
  size<-dim(X)
  n<-size[1];p<-size[2]
  
  u<-matrix(1,d+1,1)
  PotTerms<-matrix(seq(0,d),d+1,1)
  Range<-matrix(seq(0,d),d+1,1)
  
  if(p>=2){
    for (i in 2:p) {
    r1<-dim(PotTerms)[1];u1<-matrix(1,r1,1)
    ExpandedRange<-kronecker(u1,Range)
    
    PotTerms<-kronecker(PotTerms,u)
    
    PotTerms<-cbind(PotTerms,ExpandedRange)
    
  }
  }
  
  
  ActualTerms<-as.matrix(PotTerms[rowSums(PotTerms)<= d,])
  
  NumTerms<-dim(ActualTerms)[1]
  
  
  RX<-numeric()
  u<-matrix(1,p,n)
  
  for (i in 1:NumTerms) {
    
    exponent<-ActualTerms[i,]
    exponent<-t(exponent*u)
    NewTerm=apply(X^exponent,1,prod)
    
    RX<-cbind(RX,NewTerm);
    
  }
  
  return(RX)
}

ipw_fit<-function(y,D,X,taus,scale='FALSE'){
  # X are the covariates
  # D is the indicator for the respondents (treated)
  # y is the outcome
  # taus is a vector to compute quantiles for
  # Returns ATE QTE
  
  X<-as.matrix(X)
  
  y_ori=y
  if(scale=='TRUE'){y=scale(y);X=scale(X)}
  data <- data.frame(y,D,X)
  n0<-sum(D==0)
  n1<-sum(D==1)
  nt<- cbind(n0,n1)
  
  pscore.reg <- glm(D ~ X,data=data,family=binomial)
  pscore <- fitted(pscore.reg)
  
  
  weights0 <- (1-D) / (1-pscore)
  weights0 <- weights0/sum(weights0)
  weights1 <-  D / (pscore) 
  weights1 <- weights1/sum(weights1)
  weights <- weights0+weights1
  
  
  y_grid<-sort(unique(y))
  ATE <- get_ate(y,D,weights)
  cdf_qte <- get_qte(y,D,weights,taus,y_grid)
  QTE <- cdf_qte$QTE
  Q0 <- cdf_qte$Quantiles0
  Q1 <- cdf_qte$Quantiles1
  
  F0 <- cdf_qte$F0
  F1 <- cdf_qte$F1
  
  results<-list(ATE=ATE,QTE=QTE,Q0=Q0,Q1=Q1,F0 = F0,F1 = F1,weights=weights,taus=taus)
  return(results)
}

kernel_balance<-function(y,D,X,taus,nlam){
  
  Xstd <- transform.sob(X)$Xstd # standardize X to [0,1]^p
  K_X <- getGram(Xstd) # get Gram matrix using Sobolev kernel
  
  # nlam: design a grid for the tuning parameter
  lams <- exp(seq(log(1e-8), log(1), len=nlam))
  
  # compute weights for T=1
  fit1 <- ATE.ncb.SN(D, K_X, lam1s=lams)
  if (sum(fit1$warns)) cat("lambda bound warning!\n")
  
  
  # compute weights for T=0
  fit0 <- ATE.ncb.SN(1-D, K_X, lam1s=lams)
  if (sum(fit0$warns)) cat("lambda bound warning!\n")
  
  weights<-fit0$w+fit1$w
  
  y_grid<-sort(unique(y))
  ATE <- get_ate(y,D,weights)
  cdf_qte <- get_qte(y,D,weights,taus,y_grid)
  QTE <- cdf_qte$QTE
  Q0 <- cdf_qte$Quantiles0
  Q1 <- cdf_qte$Quantiles1
  
  F0 <- cdf_qte$F0
  F1 <- cdf_qte$F1
  
  results<-list(ATE=ATE,QTE=QTE,Q0=Q0,Q1=Q1,F0 = F0,F1 = F1,weights=weights,taus=taus)
  return(results)
  
}


