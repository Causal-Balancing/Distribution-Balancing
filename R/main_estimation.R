
source("R/utils.R")
source('R/ccdf_estimation.R')
source('R/weights_estimation.R')


fun_CCDF_balance_new<-function(y,D,X,taus,K=2,method_ccdf,
                           fun_obj, num_localCCDF,num_balance_grid,num_moment,
                           est_ccdf_normlize,moments_normlize,delta,
                           ...){
  X<-as.matrix(X)
  n<-length(y);p<-dim(X)[2]
  

  #Step0: randomly partition the full dataset into two equal-sized folds
  cvgroup<-make_ordered_group_balanced(n, K, D)
  
  if (num_localCCDF > 0) {
    
    # Path A: Cross-fitting with Local CCDF balancing
    # 1. Estimate CCDFs
    t1 <- proc.time()
    ccdf_results <- estimate_all_ccdfs_cv(y, D, X, cvgroup, 
                                          method_ccdf,est_ccdf_normlize, ...)
    cat(sprintf("  -> CCDF estimation done in  %.2f seconds.\n", (proc.time() - t1)[3]))
    
    # 2. Compute Weights & Counterfactuals
    t2 <- proc.time()
    counterfactual_results_list <- compute_all_counterfactuals(y, D, X, taus, cvgroup, ccdf_results,
                                                               fun_obj, num_moment, num_localCCDF,
                                                               num_balance_grid,moments_normlize,delta,...)
    cat(sprintf("  -> Weights & counterfactual distributions estimation done in %.2f seconds.\n", (proc.time() - t2)[3]))
    
    # 3. Final QTE (Calculated for Aggregated results and for each Fold)
    y_grid<-sort(unique(y))
    
    final_results_list <- lapply(counterfactual_results_list, function(results) {
      
      # --- QTE based on Aggregated Results ---
      # Method 1: Invert the Averaged CDF
      est_q0_mean <- sapply(seq_along(taus), function(i) quantile_lookup(y_grid, results$F0_hat[i,], taus[i])$quantile)
      est_q1_mean <- sapply(seq_along(taus), function(i) quantile_lookup(y_grid, results$F1_hat[i,], taus[i])$quantile)
      
      
      # --- QTE based on Fold-Specific Results ---
      # Method 1: Invert CDF for each fold
      # Loop over K folds
      est_q0_folds <- matrix(NA, nrow=K, ncol=length(taus))
      est_q1_folds <- matrix(NA, nrow=K, ncol=length(taus))
      
      for(k in 1:K) {
        est_q0_folds[k, ] <- sapply(seq_along(taus), function(i) quantile_lookup(y_grid, results$F0_hat_folds[k,i,], taus[i])$quantile)
        est_q1_folds[k, ] <- sapply(seq_along(taus), function(i) quantile_lookup(y_grid, results$F1_hat_folds[k,i,], taus[i])$quantile)
      }
      qte_cdf_folds <- est_q1_folds - est_q0_folds
      
      
      return(list(
        # Aggregated QTEs
        est_Q1=est_q1_mean,
        est_Q0= est_q0_mean,
        est_QTE_cdf_mean = est_q1_mean - est_q0_mean,
        
        # Fold-specific QTEs (Matrix: K x n_taus)
        est_QTE_cdf_folds = qte_cdf_folds))})
    
    return(list(
      path = "cross_fitting",
      cvgroup =cvgroup,
      ccdf_results = ccdf_results,
      counterfactual_results = counterfactual_results_list,
      final_results = final_results_list
    ))
    
  } else{
    # Path B: Full Sample Moment Balancing
    cat("Executing full-sample moment balancing path (num_localCCDF = 0)...\n")
    
    final_results_list <- full_sample_moment_balance(
      y, D, X, taus, fun_obj, num_moment, delta,...
    )
    
    return(list(
      path = "full_sample_moment_balancing",
      ccdf_results = NULL,
      counterfactual_results = NULL, 
      final_results = final_results_list
    ))
  }
  
}



compute_all_counterfactuals <- function(y, D, X, taus, cvgroup, ccdf_results,
                                        fun_obj, num_moment, num_localCCDF, 
                                        num_balance_grid,moments_normlize,delta, ...) {
  
  y_grid<-sort(unique(y))
  all_results <- list()

  for (moment_level in num_moment) {

    # Arrays to store CDFs for each fold
    # Dimensions: [Fold, Tau, Grid_Size]
    F0_hat_folds <- array(NA, dim = c(length(unique(cvgroup)), length(taus), length(y_grid)))
    F1_hat_folds <- array(NA, dim = c(length(unique(cvgroup)), length(taus), length(y_grid)))
    
    # list to store weights
    n_folds <- length(unique(cvgroup))
    Weights0 <- vector("list", length = n_folds)
    Weights1 <- vector("list", length = n_folds)
    
    for (k in 1:length(unique(cvgroup))) {
      
      is_target <- (cvgroup == k)
      y_target <- y[is_target]; D_target <- D[is_target]
      n_target<-length(y_target)
 
      Weights0_fold<-matrix(0,sum(1-D_target),length(taus))
      Weights1_fold<-matrix(0,sum(D_target),length(taus))
      
      Iy_target <- outer(y_target, y_grid, "<=")
      for (i in seq_along(taus)) {
        tau <- taus[i]
        # Calc Weights for D=0 and D=1
        weights0 <- get_balanced_weights(tau=tau, D_val=0, y_train=y[is_target], 
                                         D_train=D[is_target], X_train=X[is_target,], 
                                         CCDF_train=ccdf_results$CCDF0[is_target,], 
                                         y_grid=y_grid, fun_obj=fun_obj, num_moment=moment_level,
                                         num_localCCDF=num_localCCDF,
                                         num_balance_grid=num_balance_grid,
                                         moments_normlize=moments_normlize,delta,
                                          ...)
        weights1 <- get_balanced_weights(tau=tau, D_val=1, y_train=y[is_target], 
                                         D_train=D[is_target], X_train=X[is_target,], 
                                         CCDF_train=ccdf_results$CCDF1[is_target,],
                                         y_grid=y_grid, fun_obj=fun_obj, num_moment=moment_level, 
                                        num_localCCDF=num_localCCDF,  
                                         num_balance_grid=num_balance_grid,
                                        moments_normlize=moments_normlize,delta,
                                        ...)
        
        #Store Weights
        Weights0_fold[,i]<-weights0
        Weights1_fold[,i]<-weights1
        
        # Estimate Counterfactual CDFs
        
        F0_hat_folds[k, i, ] <- colSums(diag(weights0) %*% Iy_target[D[is_target] == 0,], na.rm=TRUE)
        F1_hat_folds[k, i, ] <- colSums(diag(weights1) %*% Iy_target[D[is_target] == 1,], na.rm=TRUE)
        
      }
      
      Weights0[[k]] <- Weights0_fold
      Weights1[[k]] <- Weights1_fold      
    }
    
    F0_hat_mean <- apply(F0_hat_folds, c(2, 3), mean, na.rm = TRUE)
    F1_hat_mean <- apply(F1_hat_folds, c(2, 3), mean, na.rm = TRUE)
    
    
    all_results[[paste0("moment_", moment_level)]] <- list(
      # --- Weights ---
      Weights0 = Weights0,
      Weights1 = Weights1,
      
      # --- Aggregated Results ---
      F0_hat = F0_hat_mean, 
      F1_hat = F1_hat_mean, 
      
      # --- Fold-specific Results ---
      F0_hat_folds = F0_hat_folds,     # Array: [K, n_taus, n_grid]
      F1_hat_folds = F1_hat_folds      # Array: [K, n_taus, n_grid]
    )
    
  }
 
  
  return(all_results)
}



full_sample_moment_balance <- function(y, D, X, taus, fun_obj, num_moment,delta, ...) {
  
  n <- length(y)
  all_results <- list()
  
  for (moment_level in num_moment) {
    
    cat(sprintf(" -> Processing for num_moment = %d...\n", moment_level))
    
    if (moment_level == 0) {
      q0 <- quantile(y[D == 0], taus, na.rm = TRUE)
      q1 <- quantile(y[D == 1], taus, na.rm = TRUE)
      
      all_results[[paste0("moment_", moment_level)]] <- list(est_QTE1 = q1 - q0)
      next 
    }
    
    # --- 1. 构建平衡特征 (在全样本上) ---
    moment_features <- switch(as.character(moment_level),
                              "1" = X,
                              "2" = cbind(X, X^2),
                              stop("Invalid num_moment value."))
    
    target_means <- colMeans(moment_features, na.rm = TRUE)
    
    
    # --- D=0 ---
    control_idx <- which(D == 0)
    bal_features_0 <- moment_features[control_idx, , drop = FALSE]
    base_weights_0 <- rep(1 / length(control_idx), length(control_idx))
    
    fit_w0 <- weight_w(bal_features_0, target_means, base_weights_0, obj = fun_obj,delta)
    weights0 <- fit_w0$w
    
    # --- D=1 ---
    treated_idx <- which(D == 1)
    bal_features_1 <- moment_features[treated_idx, , drop = FALSE]
    base_weights_1 <- rep(1 / length(treated_idx), length(treated_idx))
    
    fit_w1 <- weight_w(bal_features_1, target_means, base_weights_1, obj = fun_obj,delta)
    weights1 <- fit_w1$w
    
    weights<-numeric()
    weights[D==0] <-weights0
    weights[D==1] <-weights1
    
    y_grid<-sort(unique(y))
    ATE <- get_ate(y,D,weights)
    cdf_qte <- get_qte(y,D,weights,taus,y_grid)
    QTE <- cdf_qte$QTE
    Q0 <- cdf_qte$Quantiles0
    Q1 <- cdf_qte$Quantiles1
    F0 <- cdf_qte$F0
    F1 <- cdf_qte$F1
    
  
    all_results[[paste0("moment_", moment_level)]] <- list(ATE=ATE,QTE=QTE,Q0=Q0,Q1=Q1,F0 = F0,F1 = F1,weights=weights)
  }
  
  return(all_results)
}


fun_CCDF_balance_old<-function(y,D,X,taus,method_ccdf,
                               fun_obj, num_localCCDF,num_balance_grid,num_moment,
                               est_ccdf_normlize,moments_normlize,
                               ...){
  X<-as.matrix(X)
  n<-length(y);p<-dim(X)[2]
  
  
  if (num_localCCDF > 0) {
    
    # Path A: Cross-fitting with Local CCDF balancing
    # 1. Estimate CCDFs
    t1 <- proc.time()
    ccdf_results <- estimate_all_ccdfs(y, D, X, taus,method_ccdf,est_ccdf_normlize, ...)
    cat(sprintf("  -> CCDF estimation done in  %.2f seconds.\n", (proc.time() - t1)[3]))
    
    
    # 2. Compute Weights & Counterfactuals
    t2 <- proc.time()
    counterfactual_results_list <- compute_all_counterfactuals(y, D, X, taus, cvgroup=rep(1,n), ccdf_results,
                                                               fun_obj, num_moment, num_localCCDF,
                                                               num_balance_grid,moments_normlize,delta, ...)
    cat(sprintf("  -> Weights & counterfactual distributions estimation done in %.2f seconds.\n", (proc.time() - t2)[3]))
    
    # 3. Final QTE (Calculated for Aggregated results and for each Fold)
    y_grid<-sort(unique(y))
    
    final_results_list <- lapply(counterfactual_results_list, function(results) {
      
      est_q0_mean <- sapply(seq_along(taus), function(i) quantile_lookup(y_grid, results$F0_hat[i,], taus[i])$quantile)
      est_q1_mean <- sapply(seq_along(taus), function(i) quantile_lookup(y_grid, results$F1_hat[i,], taus[i])$quantile)
      
      
      return(list( est_Q1= est_q1_mean, est_Q0= est_q0_mean, est_QTE_cdf_mean = est_q1_mean - est_q0_mean))})
    
    return(list(
      path = "distribution_balancing",
      ccdf_results = ccdf_results,
      #marginal_CDF = marginal_CDF,
      counterfactual_results = counterfactual_results_list,
      final_results = final_results_list
    ))
    
  } else{
    # Path B: Full Sample Moment Balancing
    cat("Executing full-sample moment balancing path (num_localCCDF = 0)...\n")
    
    final_results_list <- full_sample_moment_balance(
      y, D, X, taus, 
      fun_obj, num_moment, geq0, sum1, ...
    )
    
    return(list(
      path = "full_sample_moment_balancing",
      ccdf_results = NULL,
      counterfactual_results = NULL, 
      final_results = final_results_list
    ))
  }
  
}





