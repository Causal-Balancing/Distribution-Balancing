
data<-read.csv("Z:\\User\\Documents\\Balance_simulations\\empirical_birth weights\\cate-birthdata\\cate-birthdata-cleaned-1stkid-white-detailed.csv")
dataset<-data[data$wt_gain!=0&data$year==2002,]

D<-dataset$smoke;y<-dataset$bweight;
sum(D);sum(D)/length(D)


#Scenorio 1
X<-cbind(dataset$mage,dataset$wt_gain)
emp_data_gen01<-function(D,X){
  
  n<-dim(X)[1]
  
  y0<-2900-4*X[,1]+X[,1]^2+5*X[,2]+0.1*X[,1]*X[,2]+1.2*X[,1]*rt(n,5)
  y1<-2600-2*X[,1]+X[,1]^2+6*X[,2]+0.2*X[,1]*X[,2]+1.5*X[,1]*rt(n,5)
  
  plot(ecdf(y0))
  lines(ecdf(y1),col='red')
  
  taus<-seq(0.05,0.95,0.01)
  QTE<-quantile(y1,taus)-quantile(y0,taus)
  plot(taus,QTE,type='l',xlab = 'Quantiles',ylab = 'True QTE')
  
  y<-D*y1+(1-D)*y0
  
  return(list(y=y,D=D,X=X,y0=y0,y1=y1))
}

set.seed(20260219)
data_1<-emp_data_gen01(D,X)
data_case1<-cbind(data_1$y,data_1$D,data_1$X,data_1$y0,data_1$y1)
write.csv(data_case1,'Z:\\User\\Documents\\Balance_simulations\\sim_empirical_birthweights\\data_case1.csv')


#Secnario 2
X<-cbind(dataset$mage,dataset$medu,dataset$X1st_prenatal,dataset$num_prenatal,
             dataset$male,dataset$married,dataset$fagemiss,dataset$diabetes,dataset$hyperpr,dataset$amnio,
             dataset$ultra,dataset$terms,dataset$drink)
    
emp_data_gen02<-function(D,X){
  
  n<-dim(X)[1]
  
  y0<-2857-0.2*X[,1]^2+14*X[,2]+27*X[,4]-19.33*X[,1]*X[,9]+1.2*X[,1]*rt(n,5)
  y0<-2392-0.3*X[,1]^2+32*X[,2]+34*X[,4]-6.61*X[,1]*X[,9]+1.5*X[,1]*rt(n,5)
  
  #plot(ecdf(y0))
  #lines(ecdf(y1),col='red')
  
  #taus<-c(0.25,0.5,0.75)
  taus<-seq(0.05,0.95,0.01)
  QTE<-quantile(y1,taus)-quantile(y0,taus)
  plot(taus,QTE,type='l',xlab = 'Quantiles',ylab = 'True QTE')
  
  y<-D*y1+(1-D)*y0
  
  return(list(y=y,D=D,X=X,y0=y0,y1=y1))
}

set.seed(20260219)
data_2<-emp_data_gen02(D,X)
data_case2<-cbind(data_2$y,data_2$D,data_2$X,data_2$y0,data_2$y1)
write.csv(data_case2,'Z:\\User\\Documents\\Balance_simulations\\sim_empirical_birthweights\\data_case2.csv')



