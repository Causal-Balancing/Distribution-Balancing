data1<-function(n){
  X<-matrix(runif(2*n,0,2),n,2)
  
  logit_p<--4+1.5*exp(X[,1]/2)+0.5*(X[,1])^2+sin(X[,1]/2)+0.4*X[,2]
  
  
  prob<-exp(logit_p)/(1+exp(logit_p))
  mean(prob);min(prob);max(prob);sum(prob<=0.1)/n;sum(prob>=0.9)/n
  
  t<-rbinom(n,1,prob)
  
  u0<-rt(n,4)
  u1<-rt(n,4)
  
  error_0<-0.3*u0*X[,2]
  error_1<-0.5*u1*X[,2]

  
  y_signal_0 <- 1.1 * X[,1]+ 0.6 * exp(X[,2])
  y_signal_1 <- 1.2 * X[,1]  +  exp(X[,2])
  #plot(density(y_signal_1),ylim=c(0,0.4));lines(density(y_signal_0),col='red')
  
  min(abs(y_signal_0))/min(abs(error_0));max(y_signal_0)/max(error_0);sd(y_signal_0)/sd(error_0)
  min(abs(y_signal_1))/mean(abs(error_1));max(y_signal_1)/max(error_1);sd(y_signal_1)/sd(error_1)
  
  
  y0<-y_signal_0+error_0
  y1<-y_signal_1+error_1
  
  y<-y1*t+y0*(1-t)
  
  return(data.frame(y=y,X=X,t=t,y0=y0,y1=y1))
}


data2<-function(n){
  X<-matrix(runif(8*n,0,2),n,8)
  
  logit_p<--4+exp(X[,1]/4)+1.5*(X[,1]/2+X[,2]/2)^2+3.5*sin(X[,2]/2)
  
  prob<-exp(logit_p)/(1+exp(logit_p))
  mean(prob);min(prob);max(prob);sum(prob<=0.1)/n;sum(prob>=0.9)/n
  
  t<-rbinom(n,1,prob)
  
  
  m0<-rt(n,4)
  m1<-rt(n,4)
  y0<-3 * exp(X[,1]) +0.5*X[,3]^2+0.7*X[,4]+0.3*m0*X[,3]
  y1<-4 * exp(X[,1]) +1.5*X[,3]^2+1.2*X[,4]+0.5*m1*X[,3]
  
  
  y<-y1*t+y0*(1-t)
  return(data.frame(y=y,X=X,t=t,y0=y0,y1=y1))
}

data3<-function(n){
  X<-matrix(runif(100*n,0,2),n,100)
  
  logit_p<--4+exp(X[,1]/4)+1.5*(X[,1]/2+X[,2]/2)^2+3.5*sin(X[,2]/2)
  
  prob<-exp(logit_p)/(1+exp(logit_p))
  mean(prob);min(prob);max(prob);sum(prob<=0.1)/n;sum(prob>=0.9)/n
  
  t<-rbinom(n,1,prob)
  
  m0<-rt(n,4)
  m1<-rt(n,4)
  y0<-3 * exp(X[,1]) +0.5*X[,3]^2+0.7*X[,4]+0.3*m0*X[,3]
  y1<-4 * exp(X[,1]) +1.5*X[,3]^2+1.2*X[,4]+0.5*m1*X[,3]
  
  
  y<-y1*t+y0*(1-t)
  return(data.frame(y=y,X=X,t=t,y0=y0,y1=y1))
}

data_ks_1<-function(n){
  Ux<-runif(n)
  Uy0<-runif(n)
  Uy1<-runif(n)
  Ud<-runif(n)
  X<-0.3+0.4*Ux
  D<-as.numeric(Ud<X)
  
  y0<-(Uy0/X)*(Uy0<=X)+Uy0*(Uy0>X)
  y1<-(Uy1/(1-X))*(Uy1<=1-X)+Uy1*(Uy1>1-X)
  
  #plot(ecdf(y0));lines(ecdf(y1),col='red')
  y<-y1*D+y0*(1-D)
  return(data.frame(y=y,X=X,t=D,y0=y0,y1=y1))
}

data_ks_2<-function(n){
  Ux<-runif(n)
  Uy0<-runif(n)
  Uy1<-runif(n)
  Ud<-runif(n)
  X<-0.3+0.4*Ux
  D<-as.numeric(Ud<X)
  
  l=0.7
  
  y0<-(Uy0<=l)*Uy0^2/l+(Uy0>l)*Uy0
  y1<-(Ux<=l)*Ux+(Ux>l)*(l+(Ux-l)^2/(1-l))
  
  #plot(ecdf(y0));lines(ecdf(y1),col='red')
  y<-y1*D+y0*(1-D)
  return(data.frame(y=y,X=X,t=D,y0=y0,y1=y1))
}


#Kallus DGP
cate_fn <- function(x){(x[,1] + x[,2]) <= 0}
makedata.sim <- function(n,d=3,ovlp=3.5){
  x <- matrix(runif(n*d),n,d)*2-1
  p <- pnorm(ovlp*((x[,1])+(x[,3]))/2)
  t <- runif(n) <= p
  cate.true <- cate_fn(x)
  f0.true <- 0.
  f1.true <- f0.true + cate.true
  eps <- (1+x[,3])*rnorm(n)
  y0 <- f0.true + eps
  y1 <- f1.true + eps
  y <- (!t)*y0 + t*y1
  return(data.frame(x=x,t,y,p,y0,y1))
}

