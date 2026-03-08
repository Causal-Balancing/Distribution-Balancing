library(quantreg)
library(quantregForest)
library(drf)
library(caret)


estimate_all_ccdfs_cv <- function(y, D, X, cvgroup, method_ccdf,est_ccdf_normlize, ...) {
  
  n <- length(y)
  y_grid_or<-sort(unique(y))
  
  if(est_ccdf_normlize=='TRUE'||est_ccdf_normlize==TRUE){
    X<-scale(X);
    y<-scale(y)
    y_grid<-sort(unique(y))
  }
  
  CCDF0 <- matrix(NA, n, length(y_grid))
  CCDF1 <- matrix(NA, n, length(y_grid))
  
  for (k in seq_len(length(unique(cvgroup)))) {
    
    is_train <- (cvgroup != k)
    is_predict  <- (cvgroup == k)
    
    # Train D=0
    train0_idx <- which(is_train & D == 0)
    
    if (length(train0_idx) > 0 && sum(is_predict) > 0) {
      
      res0 <- est_ccdf(y_train = y[train0_idx], X_train = X[train0_idx, , drop = FALSE], 
                       X_predict = X[is_predict, , drop = FALSE], y_grid = y_grid, 
                       method_ccdf = method_ccdf, est_ccdf_normlize  = est_ccdf_normlize, ...)
      
      if (!is.null(res0$CCDF)) {
        CCDF0[is_predict, ] <- res0$CCDF
      }
    }
    
    # Train D=1
    train1_idx <- which(is_train & D == 1)
    
    if (length(train1_idx) > 0 && sum(is_predict) > 0) {
      
      res1 <- est_ccdf(y_train = y[train1_idx], X_train = X[train1_idx, , drop = FALSE], 
                       X_predict = X[is_predict, , drop = FALSE], y_grid = y_grid, 
                       method_ccdf = method_ccdf, est_ccdf_normlize = est_ccdf_normlize, ...)
      
      if (!is.null(res1$CCDF)) {
        CCDF1[is_predict, ] <- res1$CCDF
      }
    }
  }
  return(list(CCDF0 = CCDF0, CCDF1 = CCDF1, y_grid = y_grid_or))
}

est_ccdf<-function(y_train,X_train,X_predict,y_grid,method_ccdf,est_ccdf_normlize,...){
  
  X_train<-as.matrix(X_train);X_predict<-as.matrix(X_predict)
  
  fit_ccdf<-switch (method_ccdf,
                    'kernel' = est_ccdf_kernel(y_train = y_train, X_train = X_train, 
                                               X_predict = X_predict, y_grid = y_grid, ...),
                     'local_linear'= est_ccdf_llr(y_train = y_train, X_train = X_train, 
                                                 X_predict = X_predict, y_grid = y_grid, ...),
                    'quantile_ref' = est_ccdf_quantile(y_train = y_train, X_train = X_train, 
                                                       X_predict = X_predict, y_grid = y_grid, ...),
                    'quantregforest' = est_ccdf_quantregforest(y_train = y_train, 
                                                               X_train = X_train, 
                                                               X_predict = X_predict,
                                                               y_grid = y_grid, ...),
                    'distforest' = est_ccdf_distforest(y_train = y_train,X_train = X_train,
                                                       X_predict = X_predict, y_grid = y_grid, ...),
                    stop(paste("Unknown method_ccdf:", method_ccdf))
  )
  
  
  return(list(CCDF=fit_ccdf$CCDF,y_grid=y_grid))
  
}

# --- Implementation of specific methods ---


est_ccdf_kernel<-function(y_train,X_train,X_predict,y_grid,kernel_type,kernel_order,
                          bandwidth_c,opt='boundary',bounds = list(c(0, 2), c(0, 2)),...){
  
  p<-dim(X_train)[2]
  #sigma<-apply(X_train,2,sd)
  #h=bandwidth_c*sigma*length(y_train)^(-1/(p+3))
  h=bandwidth_c*length(y_train)^(-1/(p+3))
  
  # --- Boundary Reflection  ---
  if (opt == 'reflection') {
  
    X_aug <- X_train
    y_aug <- y_train
    
    for (j in 1:p) {
      a <- bounds[[j]][1]
      b <- bounds[[j]][2]
      
      X_left <- X_aug
      X_left[, j] <- 2 * a - X_aug[, j]
      
      X_right <- X_aug
      X_right[, j] <- 2 * b - X_aug[, j]
      
      X_aug <- rbind(X_aug, X_left, X_right)
      y_aug <- c(y_aug, y_aug, y_aug)
    }

    X_train_final <- X_aug
    y_train_final <- y_aug
  } else {
    X_train_final <- X_train
    y_train_final <- y_train
  }
  
      
  #Step1: Kernel matriX_train
  # Use ker_h from utils.R                           
  # [nrows(X_train_final) by nrows(X_predict)]
  K_matrix <- if (missing(kernel_order)) {
    ker_h(X_train_final, X_predict, h = h, s = kernel_type)
  } else {
    Kh(X_train_final, X_predict, r = kernel_order, s = kernel_type, h = h, opt = opt)
  } 
  
  # Step 2: Indicator matrix 
  Iy <- outer(y_train_final, y_grid, "<=")    #[length(y_train) by length(y_grid)]
  
  # Step 3: Conditional CDF
  cmK <- colSums(K_matrix)       
  cmK[cmK == 0] <- .Machine$double.eps
  
  # 
  W_CCDF <- t(sweep(K_matrix, 2, cmK, "/"))   # [nrows(X_predict) by nrows(X_train_final)]
  W_CCDF[is.nan(W_CCDF)] <- 0.00001
  
  # CCDF
  CCDF <- W_CCDF %*% drop(Iy)          # [nrows(X_predict) by  length(y_grid)]        
  
  rm(K_matrix, W_CCDF)
  
  return(list(CCDF=CCDF))
  
}

est_ccdf_llr <- function(y_train, X_train, X_predict, y_grid, bandwidth_c, kernel_type = 1) {
  X_train <- as.matrix(X_train)
  X_predict <- as.matrix(X_predict)
  n <- nrow(X_train)
  nx <- nrow(X_predict)
  p <- ncol(X_train)
  
  h <- bandwidth_c *length(y_train)^(-1/(p+3))

  Iy <- outer(y_train, y_grid, "<=")  # n x J
  
  CCDF <- matrix(0, nx, length(y_grid))

  for (i in 1:nx) {
    x_i <- X_predict[i, ]

    u <- sweep(X_train, 2, x_i, "-")
    u_scaled <- sweep(u, 2, h, "/")
    
    k_vals <- 0.75 * (1 - u_scaled^2) * (abs(u_scaled) <= 1)
    weights <- apply(k_vals, 1, prod)
    
    Z <- cbind(1, u)
    idx <- which(weights > 0)
    if(length(idx) < (p + 1)) {

      CCDF[i, ] <- colMeans(Iy)
      next
    }
    
    Z_sub <- Z[idx, , drop = FALSE]
    W_sub <- weights[idx]
    Iy_sub <- Iy[idx, , drop = FALSE]
    
    ZWZ_inv <- solve(t(Z_sub * W_sub) %*% Z_sub + diag(1e-6, p+1)) 
    first_row_hat <- (ZWZ_inv %*% t(Z_sub * W_sub))[1, ]
    
    CCDF[i, ] <- first_row_hat %*% Iy_sub
  }

  CCDF[CCDF < 0] <- 0
  CCDF[CCDF > 1] <- 1
  
  return(list(CCDF = CCDF))
}

library(caret)# epsilon=0.02;m=500
est_ccdf_quantile<-function(y_train,X_train,X_predict,y_grid=y_grid,num_quantreg_grid,epsilon,...){
  
  if (!requireNamespace("quantreg", quietly = TRUE)) stop("Package 'quantreg' is required.")
  
  #Step1: Quantile regression
  
  tauseq<-seq(epsilon,1-epsilon,length=num_quantreg_grid+1)

  formula <- y_train ~ X_train 
  pred_matriX_train <- cbind(1,X_predict)

  
  bad_cols <- nearZeroVar(X_train, saveMetrics = TRUE)
  #print(bad_cols[bad_cols$nzv, ])
  X_train <- X_train[, !bad_cols$nzv]
  pred_matriX_train<-pred_matriX_train[, !bad_cols$nzv]
               
  fit <- rq(formula, tau = tauseq, data = data.frame(y_train = y_train, X_train = X_train))
  
  hatbeta <- t(fit$coefficients) # (num_quantreg_grid+1) by p or (p+1)
  
  #Step2: Estimate conditional CDF
  
  delta<-diff(tauseq)     #1 by num_quantreg_grid 
  
  pred_quantiles <- pred_matriX_train %*% t(hatbeta[-1,])   ## nrows(X_predict) by  m
  
  
  ccdf_calc <- function(pred_quantile_row, y_grid_vec, delta_vec) {
    # pred_q_row: every co in pred_quantiles, num_quantreg_grid x 1
    Iy <- outer(c(pred_quantile_row), y_grid_vec, "<=")   # num_quantreg_grid by length(y_grid)
    return(delta_vec%*% Iy)         # 1 by length(y_grid)
  }
  
  CCDF <- t(apply(pred_quantiles, 1, ccdf_calc, y_grid_vec = y_grid, delta_vec = delta))+epsilon  # nrows(X_predict) by  length(y_grid)
  
  return(list(CCDF=CCDF,hatbeta=hatbeta))
  
}


est_ccdf_quantregforest<-function(y_train,X_train,X_predict,y_grid,...){
  
  if (!requireNamespace("quantregForest", quietly = TRUE)) stop("Package 'quantregForest' is required.")
  
  fit <-quantregForest(X_train,y_train,nodesize=1, ntree=1000)
  CCDF<-predict(fit,X_predict,what=function(z){colMeans(outer(c(z),c(y_grid),'<='))}) # nrows(X_predict) by  length(y_grid)
  
  return(list(CCDF=CCDF))
  
}

est_ccdf_distforest<-function(y_train,X_train,X_predict,y_grid,min.node.size = 10,...){
  if (!requireNamespace("drf", quietly = TRUE)) stop("Package 'drf' is required.")
  
  fit <-drf(X_train,y_train,min.node.size = 10,num.trees=2000,)
  weights <- get_sample_weights(fit, newdata = X_predict)    # nrows(X_predict) by  nrows(X_train)
  
  estimate_cdf <- function(y_grid_vec, weights_matriX_train, y_train_vec) {
    
    train_y_indices <- findInterval(y_train_vec, y_grid_vec)
    
    
    CCDF <- t(apply(weights_matriX_train, 1, function(w) {
      
      weight_sums_by_interval <- tapply(w, train_y_indices, sum)
      
      full_weight_sums <- numeric(length(y_grid))
      
      valid_indices <- as.numeric(names(weight_sums_by_interval))
      
      zero_index_pos <- which(valid_indices == 0)
      if (length(zero_index_pos) > 0) {
        valid_indices <- valid_indices[-zero_index_pos]
        weight_sums_by_interval <- weight_sums_by_interval[-zero_index_pos]
      }
      
      full_weight_sums[valid_indices] <- weight_sums_by_interval
      
      cumsum(full_weight_sums)
    }))
    
    return(CCDF = CCDF)
  }
  
  CCDF <- estimate_cdf(y_grid_vec=y_grid, weights_matriX_train=weights, y_train_vec=y_train)
  
  return(list(CCDF=CCDF))
  
}

#' Full sample / Doubly Robust (DR) estimation function
estimate_all_ccdfs <- function(y, D, X, taus, method_ccdf,est_ccdf_normlize, ...) {
  
  n <- length(y)
  X<-as.matrix(X)
  y_grid<-sort(unique(y))
  y_or<- y
  y_grid_or<-y_grid
  
  data<-data.frame(y,D,X)
  if(est_ccdf_normlize==TRUE){
    X<-scale(X);
    
    mean_y0<-mean(y[D==0]);sd_y0<-sd(y[D==0])
    mean_y1<-mean(y[D==1]);sd_y1<-sd(y[D==1])
    
    y<-scale(y)
    
    y_grid<-sort(unique(y))
  }
  
  # D=0
  train0_idx <- which(D == 0)
  
  res0 <- est_ccdf(y[train0_idx], X[train0_idx, , drop = FALSE], X, y_grid, 
                   method_ccdf, est_ccdf_normlize, ... )
  CCDF0 <- res0$CCDF
  
  # D=1
  train1_idx <- which( D == 1)
  
  res1 <- est_ccdf(y[train1_idx], X[train1_idx, , drop = FALSE], X, y_grid, 
                  method_ccdf, est_ccdf_normlize, ...)
  CCDF1 <- res1$CCDF
  
  # DR Correction (Doubly Robust)
  pscore.reg <- glm(D ~ as.matrix(X),data=data,family=binomial)
  ipw <- 1/fitted(pscore.reg)
  
  Iy<-matrix(outer(y_or,y_grid_or,"<="),n,length(y_grid_or))
  ipw0_matrix<-matrix(ipw*(1-D)/sum(ipw*(1-D)),n,length(y_grid_or),byrow=FALSE)
  hat_cdf0_dr<-colMeans(CCDF0+(Iy-CCDF0)*ipw0_matrix)
  
  ipw1_matrix<-matrix(ipw*D/sum(ipw*D),n,length(y_grid_or),byrow=FALSE)
  hat_cdf1_dr<-colMeans(CCDF1+(Iy-CCDF1)*ipw1_matrix)
  
  # Lookup Quantiles
  if(est_ccdf_normlize==TRUE){
    Q1<-quantile_lookup(y_grid_or,colMeans(CCDF1), taus)$quantile
    Q0<-quantile_lookup(y_grid_or,colMeans(CCDF0), taus)$quantile
    
    Q1_DR<-quantile_lookup(y_grid_or,hat_cdf1_dr, taus)$quantile
    Q0_DR<-quantile_lookup(y_grid_or,hat_cdf0_dr, taus)$quantile
  }else{
    Q1<-quantile_lookup(y_grid,colMeans(CCDF1), taus)$quantile
    Q0<-quantile_lookup(y_grid,colMeans(CCDF0), taus)$quantile
    
    Q1_DR<-quantile_lookup(y_grid,hat_cdf1_dr, taus)$quantile
    Q0_DR<-quantile_lookup(y_grid,hat_cdf0_dr, taus)$quantile
  }
  
  return(list(CCDF0 = CCDF0, CCDF1 = CCDF1, QTE = Q1-Q0,Q1=Q1,Q0=Q0,
              CDF0_DR = hat_cdf0_dr, CDF1_DR = hat_cdf1_dr,QTE_DR =Q1_DR- Q0_DR,
              y_grid = y_grid))
}


