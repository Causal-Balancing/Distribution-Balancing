DH_test<-function(y,D,X,d,B,alpha){
  
  X<-as.matrix(X)
  size<-dim(X)
  n<-size[1];p<-size[2]
  
  RX<-powers_of_x(X,d)
  pscore.reg <- glm(D ~ RX,family=binomial)
  pscore <- fitted(pscore.reg)
  
  weight00<-(1-D)/(1-pscore)
  weight0<-weight00/sum(weight00)
  
  weight11<-D/pscore
  weight1<-weight11/sum(weight11)
  
  q<-unique(sort(y));m<-length(q)
  I0<-matrix(1,n,m);I1<-matrix(1,m,n)
  Ym=c(y)*I0;Qm=t(q*I1)
  Iy<-(Ym<=Qm) # n*m
  
  F0hat<-c(weight0)%*%Iy
  F1hat<-c(weight1)%*%Iy
  
  KShat<-sqrt(n)*max(F1hat-F0hat)
  
  F0xhat<-t(weight00*Iy)%*%RX%*%solve(t(RX)%*%RX)%*%t(RX) #m*n
  
  F1xhat<-t(weight1*Iy)%*%RX%*%solve(t(RX)%*%RX)%*%t(RX)
 
  
  S<-numeric()
  for (b in 1:B) {
    U<-rnorm(n,0,1)
    Psi0<-t(U*weight00)%*%Iy-F0hat*sum(U)+t(U*(1-weight00))%*%t(F0xhat)
    Psi0<-Psi0/sqrt(n)
    
    
    Psi1<-t(U*weight11)%*%Iy-F0hat*sum(U)+t(U*(1-weight11))%*%t(F0xhat)
    Psi1<-Psi1/sqrt(n)
    
    S[b]<-max(Psi1-Psi0)
    
  }
  
  chat<-numeric()
  rejection<-numeric()
  for (a in 1:length(alpha)) {
    c<-quantile(S,1-alpha[a])
    chat[a]<-c
    rejection[a]<-(KShat>c)
  }

  #rejection<-(KShat>chat)
 
  return(list(KShat=KShat,chat=chat,alpha=alpha,rejection=rejection))
}

SLE_test<-function(y,D,X,d,B,alpha){
  
  dataset<-cbind(y,D,X)
  
  n<-length(y)
  y_grid<-sort(unique(y))
  
  fit<-SLE(y,D,X,d=2,taus)
  F1_hat<-fit$F1
  F0_hat<-fit$F0
  
  KS_hat<-sqrt(n)*max(F1_hat-F0_hat)
  
  KS_hat_boot<-numeric()
  
  for (b in 1:B) {
    index<-sample(1:n,n,replace = TRUE)
    
    data_boot<-dataset[index,]
    
    y_boot<-data_boot[,1];
    D_boot<-data_boot[,2]
    X_boot<-data_boot[,-c(1,2)]
    
    fit_SLE_boot<-SLE(y_boot,D_boot,X_boot,d=2,taus,y_grid = y_grid)
    
    F1_hat_boot<-fit_SLE_boot$F1
    F0_hat_boot<-fit_SLE_boot$F0
    
    KS_hat_boot[b]<-sqrt(n)*max(F1_hat_boot-F0_hat_boot-(F1_hat-F0_hat))
  }
  
  pvalue<- (sum(KS_hat_boot > KS_hat)+1)/(B+1)
  
  rejection<-numeric()
  for (a in 1:length(alpha)) {
    
    rejection[a]<-as.numeric(pvalue<alpha[a])
  }
  
  
  return(list(KS_hat=KS_hat,pvalue=pvalue,rejection=rejection))
  
}
