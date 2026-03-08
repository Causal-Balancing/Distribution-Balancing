
#' Generate cross-validation groups
make_cvgroup = function(n, K, right = TRUE) {
  split     = runif(n)
  return(as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)), include.lowest = TRUE, right =TRUE)))
}

#' Generate balanced cross-validation groups (stratified by Treatment)
make_cvgroup_balanced = function(n, K, D) {
  cvgroup = numeric(n)
  cvgroup[D==1] = make_cvgroup(sum(D==1), K, right = TRUE)
  cvgroup[D==0] = make_cvgroup(sum(D==0), K, right = TRUE)
  return(cvgroup)
}

make_ordered_group = function(n, K = 2) {
  sub_size = floor(n / K)
  return(c(rep(1, sub_size), rep(2, n - sub_size)))
}

make_ordered_group_balanced = function(n, K, D) {
  split_group = numeric(n)
  split_group[D == 1] = make_ordered_group(sum(D == 1), K = K)
  split_group[D == 0] = make_ordered_group(sum(D == 0), K = K)
  return(split_group)
}

get_ate<-function(y,D,weights){
  weights0<-weights[D==0]
  weights1<-weights[D==1]
  ate <- getWeightedMean(y[D==1],weights1)-getWeightedMean(y[D==0],weights0)
  return(ate)
}

get_cdf<-function(y,weights,y_grid){
  weights<-weights/sum(weights)
  Iy<- outer(y, y_grid, "<=")
  CDF <- colSums(diag(weights) %*% Iy)
  return(CDF = CDF)
}

get_qte<-function(y,D,weights,taus,y_grid){
  weights0<-weights[D==0]
  weights1<-weights[D==1]
  
  F0<-get_cdf(y[D==0],weights0,y_grid)
  F1<-get_cdf(y[D==1],weights1,y_grid)
  
  Quantiles0<-sapply(seq_along(taus), function(i) quantile_lookup(y_grid, F0, taus[i])$quantile)
  Quantiles1<-sapply(seq_along(taus), function(i) quantile_lookup(y_grid, F1, taus[i])$quantile)
  
  QTE <- Quantiles1 - Quantiles0
  names(QTE)<-c(taus)
  return(list(F0 = F0,F1 = F1, Quantiles1=Quantiles1,Quantiles0=Quantiles0,QTE = QTE))
}




#' Inverse lookup of quantiles from CDF
quantile_lookup <- function(y_grid, marginal_cdf, taus) {
  phi<-marginal_cdf/max(marginal_cdf)
  marginal_cdf<-as.numeric(lapply(y_grid,function(x) max(0,max(phi[y_grid<=x]))))
  
  quantile<-sapply(taus, function(t) min(y_grid[marginal_cdf>=t]) )
  
  return(list(marginal_cdf=marginal_cdf,quantile=quantile))
}

#' Generalized kernel function
ker<-function(u,s){
  u<-as.matrix(u)
  
  poch<-function(a,n){
    gamma(a+n)/gamma(a)
  }
  
  if(is.finite(s) & s>=0 & s==floor(s)){
    K=poch(1/2,s+1)/prod(1:s)*(1-u^2)^s*(abs(u)<=1);
  }else if(is.infinite(s)){
    K = (2*pi)^(-1/2)*exp(-u^2/2);
  }else {
    stop('s must be a nonnegative integer or s = inf for Gaussian kernel.')
  }
  return(K)
}

#' Multivariate kernel matrix computation
ker_h<-function(X,x,h,s){
  X<-as.matrix(X)
  n<-dim(X)[1];p<-dim(X)[2]
  if(is.null(dim(x))){
    px<-length(x);nx=1
    x<-matrix(x,1,px)
  }else{
    nx<-dim(x)[1];px<-dim(x)[2]
  }
  
  if(p!=px){
    stop("X and x must be of the same dimension.")
  }
  
  dist<-rep(NULL,p)
  cont<-rep(NULL,p)
  for (j in 1:p) {
    dist[j]<-all(X[,j]==0 | X[,j]==1) & any(X[,j]!=0) & any(X[,j]!=1)
    cont[j]<-any(X[,j]!=0 & X[,j]!=1) | all(X[,j]==0) | all(X[,j]==1)
  }
  
  if(sum(dist)==0){
    Xd=NULL;xd=NULL;Xc=X;xc=x
  }else if(sum(cont)==0){
    Xc=NULL;xc=NULL;Xd=X;xd=x
  }else{
    Xd = matrix(X[,dist],n,sum(dist))
    xd = matrix(x[,dist],nx,sum(dist))
    Xc = matrix(X[,cont],n,sum(cont))
    xc = matrix(x[,cont],nx,sum(cont))
  }
  
  pc<- dim(Xc)[2];pd<-dim(Xd)[2]
  
  if(length(h)==1){h = h*rep(1,pc)}
  
  Kernel<-matrix(0,n,nx)
  
  if(sum(dist)==0){
    for (i in 1:n) {
      ac<-Xc[i,]
      Ac<-matrix(ac,nx,pc,byrow = TRUE)
      scaled_diff <- sweep(Ac - xc, 2, h, "/")  
      k_h <- sweep(ker(scaled_diff, s), 2, h, "/")
      Kernel[i,]<-apply(k_h, 1, prod)
    }
  }else if(sum(cont)==0){
    for (i in 1:n) {
      ad <- Xd[i, ]
      Ad <- matrix(ad, nx, pd, byrow = TRUE)
      Kernel[i, ] <- apply(Ad == xd, 1, prod)
    }
  }else{
    for (i in 1:n){
      ac <- Xc[i, ]
      Ac <- matrix(ac, nx, pc, byrow = TRUE)
      u <- sweep(Ac - xc, 2, h, "/")
      k_val <- sweep(ker(u, s), 2, h, "/")
      k_cont <- apply(k_val, 1, prod)

      ad <- Xd[i, ]
      Ad <- matrix(ad, nx, pd, byrow = TRUE)
      k_dist <- apply(Ad == xd, 1, prod)
      
      Kernel[i, ] <- k_cont * k_dist
    }
  }
  
  return(Kernel)
}#gen_kernel_mul(U-X,h,type)


# Higher-order smooth polynomial kernels
#
# K(u,r,s) returns an (n x p) rth-order, s-smooth polynomial kernel from an (n x k) matrix u. 
# u: (n x k) matrix 
# The order can be any positive even number.
# smoothness index s can be any non-negative integer: 0 = uniform, 1 =
# Epanechnikov, 2 = biweight, 3 = triweight, and inf = Gaussian kernel. See
# Wand and Schucany (1990) and Hansen (2005) for more details.

fun_K<-function(u,r,s){
  poch<-function(a,n){
    gamma(a+n)/gamma(a)
  }
  if(is.null(dim(u))){
    nu<-1;pu<-length(u)
  }else{
    nu<-dim(u)[1];pu<-dim(u)[2]
  }
  
  if(r>0 & r/2==floor(r/2)){
    r<-r/2
    c<-array(0,dim=c(nu,pu,r))
    if(is.finite(s) & s>=0 & s==floor(s)){
      for (i in 1:r) {
        c[,,i]<-(-1)^(i-1)*poch(1/2+s+r,i-1)*u^(2*(i-1))/(factorial(i-1)*factorial(r-i)*poch(3/2,i-1))
      }
      B = poch(3/2,r-1)*poch(3/2+s,r-1)/poch(s+1,r-1)*apply(c,c(1,2),sum);
      M = poch(1/2,s+1)/prod(1:s)*(1-u^2)^s*(abs(u)<=1);
      out = B*M;
    }else if(is.infinite(s)){
      for (i in 1:r) {
        c[,,i] = (-1)^(i-1)*2^(i-2*r)*factorial(2*r)*u^(2*(i-1))/(factorial(r)*factorial(2*i-1)*factorial(r-i));
      }
      phi = (2*pi)^(-1/2)*exp(-u^2/2);
      out = apply(c,c(1,2),sum)*phi;
    }else {
      stop('s must be a nonnegative integer or s = inf for Gaussian kernel.')
    }
  }else{
    stop('r must be a positive even number.')
  }
  return(out)
}

# Higher-order boundary kernels
#
# Kh(X,x,r,s,h,opt) returns an (n x nx) higher-order boundary kernel from an
# (n x p) data matrix X and an (nx x px) matrix of evaluation points x. 
# The boundary region is determined by the bandwidth vector h for continuous
# variables in X. In addition, Kh(X,x,r,s,h,opt) handles the discrete
# variables using frequency-based method as in Li and Racine (2007).
#
# Kernel: (n x nx) output
# X: (n x p) covariates; can be continuous, discrete, or mixed
# x: (nx x px) evaluation points
# r: Kernel order; must be an even number
# s: Kernel smoothness: 0=uniform; 1=Epanechnikov; 2=biweight; 3=triweight;
#    inf=Gaussian
# h: (1 x kc) bandwidths for continuous variables in X. Specify h=[] for 
#    the rule-of-thumb bandwidths
# opt: [] or 'boundary' for boundary kernel; 'standard' for stnadard kernel

Kh<-function(X,x,r,s,h,opt){
  
  X<-as.matrix(X)
  n<-dim(X)[1];p<-dim(X)[2]
  if(is.null(dim(x))){
    px<-length(x)
    nx=1
    x<-matrix(x,1,px)
  }else{
    nx<-dim(x)[1];px<-dim(x)[2]
  }
  
  if(p!=px){
    stop("X and x must be of the same dimension.")
  }
  
  dist<-rep(NULL,p)
  cont<-rep(NULL,p)
  for (j in 1:p) {
    dist[j]<-all(X[,j]==0 | X[,j]==1) & any(X[,j]!=0) & any(X[,j]!=1)
    cont[j]<-any(X[,j]!=0 & X[,j]!=1) | all(X[,j]==0) | all(X[,j]==1)
  }
  
  if(sum(dist)==0){
    Xd=NULL;xd=NULL;Xc=X;xc=x
  }else if(sum(cont)==0){
    Xc=NULL;xc=NULL;Xd=X;xd=x
  }else{
    Xd = as.matrix(X[,dist],n,sum(dist))
    xd = as.matrix(x[,dist],nx,sum(dist))
    Xc = as.matrix(X[,cont],n,sum(cont))
    xc = as.matrix(x[,cont],nx,sum(cont))
  }
  
  
  pc<- dim(Xc)[2];pd<-dim(Xd)[2]
  
  fun1<-function(u){
    fun_K(u,r,s)^2
  }
  fun2<-function(u){
    u^r*fun_K(u,r,s)
  }
  if(missing(h)){ #rule of thumb bandwidth
    
    R<-integrate(fun1,-Inf,Inf)$value
    kappa<-integrate(fun2,-Inf,Inf)$value
    
    fun_of<-function(n){
      factorial(n+1)/(2^((n+1)/2)*factorial((n+1)/2));
    }
    c = (pi^(pc/2)*2^(pc+r-1)*factorial(r)^2*R^pc/(r*kappa^2*(fun_of(2*r-1)+(pc-1)*fun_of(r-1)^2)))^(1/(2*r+pc));
    h = c*apply(Xc,2,sd)*n^(-1/(2*r+pc));
    
  }else if(length(h)==1){
    h = h*rep(1,pc);
  }else if(length(h)!=1 & length(h)!=pc && any(h<=0)){
    stop('h must be a positive scalar or a vector of the same dimension as continuous variables in X.');
  }
  
  Kernel<-matrix(0,n,nx) #Output
  
  divide<-function(A,a){
    t(apply(A, 1, function(z) z/a ))
  }
  
  if(opt=='boundary'){ #% use local linear/quadratic fit to construct boundary kernel
    
    Xc_min<-matrix(apply(Xc,2,min),nx,p,byrow = TRUE)
    dl = divide(Xc_min-xc,h)
    dl[dl<=-1]=-1
    
    Xc_max<-matrix(apply(Xc,2,max),nx,p,byrow = TRUE)
    dh =divide(Xc_max-xc,h)
    dh[dh>=1]=1
    
    mu0<-matrix(0,nx,p);mu1<-matrix(0,nx,p);mu2<-matrix(0,nx,p);
    for (i in 1:nx) {
      for (j in 1:p) {
        fun3<-function(t){
          fun_K(dl[i,j]+t*(dh[i,j]-dl[i,j]),2,s)*(dh[i,j]-dl[i,j])
        }
        fun4<-function(t){
          (dl[i,j]+t*(dh[i,j]-dl[i,j]))*fun_K(dl[i,j]+t*(dh[i,j]-dl[i,j]),2,s)*(dh[i,j]-dl[i,j])
        }
        
        fun5<-function(t){
          (dl[i,j]+t*(dh[i,j]-dl[i,j]))^2*fun_K(dl[i,j]+t*(dh[i,j]-dl[i,j]),2,s)*(dh[i,j]-dl[i,j])
          
        }
        
        mu0[i,j]<-integrate(fun3,0,1)$value
        mu1[i,j]<-integrate(fun4,0,1)$value
        mu2[i,j]<-integrate(fun5,0,1)$value
        
      }
      
    }
    #mu0[mu0==0]<-1e-20
    #mu1[mu1==0]<-1e-20
    #mu1[mu1==0]<-1e-20
    k = r-pc;
    if (k <= 0){ #local constant
      fun_BK = function(u) {mu0^(-1)*fun_K(u,2,s)}
    }else if(k == 1 ){ #local linear 
      fun_BK = function(u){(mu2-mu1*u)*fun_K(u,2,s)/(mu0*mu2-mu1^2)}
    }else if (k==2){ #local quadratic
      
      mu3<-matrix(0,nx,p);mu4<-matrix(0,nx,p);
      for (i in 1:nx) {
        for (j in 1:p) {
          fun6<-function(t){
            (dl[i,j]+t*(dh[i,j]-dl[i,j]))^3*fun_K(dl[i,j]+t*(dh[i,j]-dl[i,j]),2,s)*(dh[i,j]-dl[i,j])
          }
          fun7<-function(t){
            (dl[i,j]+t*(dh[i,j]-dl[i,j]))^4*fun_K(dl[i,j]+t*(dh[i,j]-dl[i,j]),2,s)*(dh[i,j]-dl[i,j])
          }
          
          mu3[i,j]<-integrate(fun6,0,1)$value
          mu4[i,j]<-integrate(fun7,0,1)$value
          
        }
      }
      
      
      
      fun_BK<-function(u){
        (mu2*mu4-mu3^2+(mu2*mu3-mu1*mu4)*u+(mu1*mu3-mu2^2)*u^2)*fun_K(u,2,s)/
          (mu0*mu2*mu4+2*mu1*mu3*mu2-mu2^3-mu1^2*mu4-mu0*mu3^2);
      }
    }else{
      stop('In the boundary kernel case, kernel order must be the smallest even number greater than the number of continuous variables.');
      
    }
    
    if(sum(dist)==0){
      for(i in 1:n){
        ac<-Xc[i,]
        Ac<-matrix(ac,nx,pc,byrow = TRUE)
        #apply(xc, 1,function(x) Xc[i,]-x )
        Uc<-divide(Ac-xc,h)
        Kernel[i,] = t(apply(divide(fun_BK(Uc),h),1,prod))
      }
    }else if(sum(cont)==0){
      for (i in 1:n){
        ad<-Xd[i,]
        Ad<-matrix(ad,nx,pd,byrow = TRUE)
        Ud<-divide(Ad-xd,h)
        Kernel[i,] = t(apply(divide(fun_BK(Ud),h),1,prod));
      }
    }else{
      for (i in 1:n){
        ac<-Xc[i,]
        Ac<-matrix(ac,nx,pc,byrow = TRUE)
        ad<-Xd[i,]
        Ad<-matrix(ad,nx,pd,byrow = TRUE)
        Uc<-divide(Ac-xc,h)
        Ud<-divide(Ad-xd,h)
        Kernel[i,] = t(apply(divide(fun_BK(Uc),h),1,prod))*t(apply(divide(fun_BK(Ud),h),1,prod))
      }
    }
  }else if(opt=='standard'){
    if(sum(dist)==0){
      for (i in 1:n){
        ac<-Xc[i,]
        Ac<-matrix(ac,nx,pc,byrow = TRUE)
        #apply(xc, 1,function(x) Xc[i,]-x )
        Uc<-divide(Ac-xc,h)
        Kernel[i,] = t(apply(divide(fun_K(Uc,r,s),h),1,prod))
      }
    }else if(sum(cont)==0){
      for (i in 1:n){
        ad<-Xc[i,]
        Ad<-matrix(ad,nx,pc,byrow = TRUE)
        #apply(xc, 1,function(x) Xc[i,]-x )
        Ud<-divide(Ad-xd,h)
        Kernel[i,] = t(apply(divide(fun_K(Ud,r,s),h),1,prod))
      }
    }else{
      for (i in 1:n){
        ac<-Xc[i,]
        Ac<-matrix(ac,nx,pc,byrow = TRUE)
        ad<-Xd[i,]
        Ad<-matrix(ad,nx,pd,byrow = TRUE)
        Uc<-divide(Ac-xc,h)
        Ud<-divide(Ad-xd,h)
        Kernel[i,] = t(apply(divide(fun_K(Uc,r,s),h),1,prod))*t(apply(divide(fun_K(Ud,r,s),h),1,prod))
        
      }
    }
    
  }else{
    stop('opt must be boundary or standard.')
  }
  
  return(Kernel)
  
}
