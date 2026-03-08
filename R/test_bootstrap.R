fun_CCDF_balance_with_test <- function(y, D, X, taus, K = 2, method_ccdf,
                                            fun_obj, num_localCCDF, num_balance_grid, num_moment,
                                            est_ccdf_normlize, moments_normlize,
                                            num_boot = 500, alpha = 0.05,delta = 0.01,
                                            bootstrap_dist = "exponential", seed = NULL,
                                            maxit = 500,
                                            ...) {
  fit <- fun_CCDF_balance_new(
    y = y, D = D, X = X, taus = taus, K = K,
    method_ccdf = method_ccdf,
    fun_obj = fun_obj,
    num_localCCDF = num_localCCDF,
    num_balance_grid = num_balance_grid,
    num_moment = num_moment,
    est_ccdf_normlize = est_ccdf_normlize,
    moments_normlize = moments_normlize,delta,
    ...
  )
  
  test_infer <- sd_test_multiplier(
    y = y, D = D, X = X, taus = taus,
    fit_obj = fit,
    num_moment = num_moment,
    num_localCCDF = num_localCCDF,
    num_balance_grid = num_balance_grid,
    fun_obj = fun_obj,
    moments_normlize = moments_normlize,
    num_boot = num_boot,
    alpha = alpha,
    delta = delta,
    bootstrap_dist = bootstrap_dist,
    seed = seed,
    maxit = maxit
  )
  
  fit$test <- test_infer
  fit
}

sd_test_multiplier <- function(y, D, X, taus, fit_obj,
                                     num_moment, num_localCCDF, num_balance_grid,
                                     fun_obj = "entropy",
                                     moments_normlize,
                                     num_boot = 500, alpha = 0.1,delta = 0.01,
                                     bootstrap_dist = "exponential",
                                     seed = NULL,
                                     maxit = 200) {
  if (!is.null(seed)) set.seed(seed)
  if (fit_obj$path != "cross_fitting") {
    stop("Multiplier bootstrap is implemented for cross-fitting path only (num_localCCDF > 0).")
  }
  if (fun_obj != "entropy") {
    stop("Current multiplier bootstrap implementation supports fun_obj='entropy' only.")
  }
  
  X <- as.matrix(X)
  cvgroup <- fit_obj$cvgroup
  ccdf_results <- fit_obj$ccdf_results
  y_grid <- ccdf_results$y_grid
  K <- length(unique(cvgroup))
  n<-length(y)
  out <- list()
  
  for (moment_level in num_moment) {
    test_statistics_boot <- numeric()
    F0_hat <- fit_obj$counterfactual_results[[paste0("moment_", moment_level)]]$F0_hat
    F1_hat <- fit_obj$counterfactual_results[[paste0("moment_", moment_level)]]$F1_hat
    test_statistics <- sqrt(n)*max(F1_hat -F0_hat)
    
    designs_precomputed <- list()
    for (i in seq_along(taus)) {
      tau <- taus[i]
      designs_precomputed[[i]] <- list()
      for (k in seq_len(K)) {
        is_target <- (cvgroup == k)
        
        d0 <- build_balance_design(tau, y, D, X, ccdf_results$CCDF0, y_grid, num_localCCDF, num_balance_grid, moment_level, moments_normlize)
        d1 <- build_balance_design(tau, y, D, X, ccdf_results$CCDF1, y_grid, num_localCCDF, num_balance_grid, moment_level, moments_normlize)
        
        designs_precomputed[[i]][[k]] <- list(d0 = d0, d1 = d1)
      }
    }
    
    for (b in 1:num_boot) {
      
      xi <- .multiplier_draw(length(y), dist = bootstrap_dist)
      
      F0_boot_folds <- array(NA_real_, dim = c(K, length(taus), length(y_grid)))
      F1_boot_folds <- array(NA_real_, dim = c(K, length(taus), length(y_grid)))
      
      tau_boot_results <- numeric(length(taus))
      
      for (i in seq_along(taus)) {
        F_fold_0 <- matrix(NA_real_, K, length(y_grid))
        F_fold_1 <- matrix(NA_real_, K, length(y_grid))
        
        for (k in seq_len(K)) {
          is_target <- (cvgroup == k)
          idx0 <- which(is_target & D == 0)
          idx1 <- which(is_target & D == 1)
          
          
          d0_all <- designs_precomputed[[i]][[k]]$d0
          d1_all <- designs_precomputed[[i]][[k]]$d1
          
          
          M0_boot <- colMeans(d0_all$bal_features[is_target,] * xi[is_target]) 
          M1_boot <- colMeans(d1_all$bal_features[is_target,] * xi[is_target]) 
          
          C0 <- d0_all$bal_features[idx0, , drop = FALSE]
          C1 <- d1_all$bal_features[idx1, , drop = FALSE]
          
          
          w0 <- weight_entropy_multipiler(C0, M0_boot, rep(1, length(idx0)), xi[idx0], delta, maxit)
          w1 <- weight_entropy_multipiler(C1, M1_boot, rep(1, length(idx1)), xi[idx1], delta, maxit)
          
          
          Iy0 <- outer(y[idx0], y_grid, "<=")
          Iy1 <- outer(y[idx1], y_grid, "<=")
          F_fold_0[k, ] <- drop(crossprod(w0, Iy0))
          F_fold_1[k, ] <- drop(crossprod(w1, Iy1))
        }
        
        
        F0_final <- colMeans(F_fold_0)
        F1_final <- colMeans(F_fold_1)
        
        test_statistics_boot_results<-sqrt(n) * max(F1_final - F0_final-(F1_hat-F0_hat))
        
      }
      test_statistics_boot[b] <- test_statistics_boot_results
    }
    
    pvalue<- (sum(test_statistics_boot > test_statistics)+1)/(num_boot+1)
    
    out[[paste0("moment_", moment_level)]] <- list(
      test_statistics=test_statistics,
      pvalue = pvalue,
      alpha = alpha,
      bootstrap_dist = bootstrap_dist
    )
  }
  
  out
}



