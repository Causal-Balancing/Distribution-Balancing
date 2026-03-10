rm(list=ls())
gc()

# ====================1. Environment Setup & Libraries===============================
packages.needed <- c('foreach', 'parallel', 'doParallel', 'doRNG', 'quantreg', 'drf', 
                     'quantregForest', 'osqp','BMisc')
lapply(packages.needed, library, character.only = TRUE)

# Load core functions (assuming current working dir is the project root)
source('R/main_estimation.R')


# ================================2. Data====================================
data<-read.csv("Applications/cate-birthdata-cleaned-1stkid-white-detailed.csv")

dataset<-data[data$wt_gain!=0&data$year==2002,]
write.csv(dataset,'dataset.csv')
D<-dataset$smoke;y<-dataset$bweight;
X<-cbind(dataset$wt_gain,dataset$mage,dataset$medu,dataset$X1st_prenatal,dataset$num_prenatal,
         dataset$male,dataset$married,dataset$fagemiss,dataset$diabetes,dataset$hyperpr,dataset$amnio,
         dataset$ultra,dataset$terms,dataset$drink)
mydata<-data.frame(y,D,X)

# ================================3. Estimation====================================

taus<-seq(0.1,0.9,0.05)

library(recipes)
if(is.null(colnames(X))) {
  colnames(X) <- paste0("x", 1:14)
}
my_data <- as.data.frame(X)

#Expantion X
rec <- recipe(~ ., data = my_data) %>%
  step_poly(all_numeric_predictors(), degree = 2, options = list(raw = TRUE)) %>%
  step_interact(terms = ~ all_predictors():all_predictors()) %>%
  step_nzv(all_predictors()) %>% 
  step_lincomb(all_predictors()) 

prep_data <- prep(rec)
final_X <- bake(prep_data, new_data = my_data)

# 'quantile_ref'
# 'quantregforest'
# 'distforest'
#  distribution balancing
fit_ccdf_balance <- fun_CCDF_balance_new(y, D, final_X, taus,K=2, method_ccdf='distforest',
                                         fun_obj='entropy', num_localCCDF=2, num_balance_grid = 100,
                                         num_moment=c(1), 
                                         moments_normlize = TRUE,
                                         est_ccdf_normlize=TRUE, 
                                         num_quantreg_grid=200,epsilon=0.05,delta =0.01)

QTE_Epy3<-fit_ccdf_balance$final_results$moment_1$est_QTE_cdf_mean
plot(taus,QTE_Epy3,type='l')

F1_hat <- fit_ccdf_balance$counterfactual_results$moment_1$F1_hat[1,]
F0_hat <- fit_ccdf_balance$counterfactual_results$moment_1$F0_hat[1,]
n <- length(y)
Test_statistics <- sqrt(n)*max(F0_hat-F1_hat)


#ATE
#library("grf")
rf_prob <- probability_forest(X = as.matrix(final_X), Y = as.factor(D), num.trees = 2000)
pi_hat <- predict(rf_prob)$predictions[, 2]
pi_hat <- pmin(pmax(pi_hat, 0.01), 0.99)
w1 <- D / pi_hat
y1_hat <- sum(w1 * y) / sum(w1)
w0 <- (1 - D) / (1 - pi_hat)
y0_hat <- sum(w0 * y) / sum(w0)

ate_hat <- y1_hat - y0_hat
  

# ================================3. Inference====================================

cvgroup <- fit_ccdf_balance$cvgroup
ccdf_results <- fit_ccdf_balance$ccdf_results
y_grid <- sort(unique(y))
K <- length(unique(cvgroup))
moment_level <- 1
qte_point <- fit_ccdf_balance$final_results[[paste0("moment_", moment_level)]]$est_QTE_cdf_mean

designs_precomputed <- list()
for (i in seq_along(taus)) {
  tau <- taus[i]
  designs_precomputed[[i]] <- list()
  for (k in seq_len(K)) {
    is_target <- (cvgroup == k)
    
    d0 <- build_balance_design(tau, y, D, final_X, ccdf_results$CCDF0, y_grid, num_localCCDF=3, num_balance_grid=100, num_moment=moment_level, moments_normlize=TRUE)
    d1 <- build_balance_design(tau, y, D, final_X, ccdf_results$CCDF1, y_grid, num_localCCDF=3, num_balance_grid=100, num_moment=moment_level, moments_normlize=TRUE)
    
    designs_precomputed[[i]][[k]] <- list(d0 = d0, d1 = d1)
  }
}

# ==============================4. Parallel =============================

detectCores(no_cores)
no_cores<-80

cl = parallel::makeCluster(no_cores, type = "SOCK")
registerDoParallel(cl)

# Set random seed stream for reproducibility
parallel::clusterSetRNGStream(cl,)


# Export Environment to Cluster
clusterEvalQ(cl, {
  library(quantreg)
  library(BMisc)
  library(quantregForest)
  library(drf)
  library(osqp)
  source('R/main_estimation.R')
  source('R/inference_bootstrap.R')
  source('sim_utils/SLE.R')
  source('sim_utils/LDML.R')
})

num_boot<-1000

start_time <- proc.time() 
result <- foreach(k = 1:num_boot, .combine = 'rbind',
                  .packages = c('quantreg', 'drf', 'quantregForest', 'osqp')) %dopar% {
                    
                    
                    xi <- .multiplier_draw(length(y), dist = 'exponential')
                    
                    
                    F_fold_0 <- matrix(NA_real_, K, length(y_grid))
                    F_fold_1 <- matrix(NA_real_, K, length(y_grid))
                    
                    for (k in seq_len(K)) {
                      is_target <- (cvgroup == k)
                      idx0 <- which(is_target & D == 0)
                      idx1 <- which(is_target & D == 1)
                      
                      
                      d0_all <- designs_precomputed[[1]][[k]]$d0
                      d1_all <- designs_precomputed[[1]][[k]]$d1
                      
                      
                      M0_boot <- colMeans(d0_all$bal_features[is_target,] * xi[is_target]) 
                      M1_boot <- colMeans(d1_all$bal_features[is_target,] * xi[is_target]) 
                      
                      C0 <- d0_all$bal_features[idx0, , drop = FALSE]
                      C1 <- d1_all$bal_features[idx1, , drop = FALSE]
                      
                      
                      w0 <- weight_entropy_multipiler(C0, M0_boot, rep(1, length(idx0)), xi[idx0], delta=0.01, maxit=500)
                      w1 <- weight_entropy_multipiler(C1, M1_boot, rep(1, length(idx1)), xi[idx1], delta=0.01, maxit=500)
                      
                      
                      Iy0 <- outer(y[idx0], y_grid, "<=")
                      Iy1 <- outer(y[idx1], y_grid, "<=")
                      F_fold_0[k, ] <- drop(crossprod(w0, Iy0))
                      F_fold_1[k, ] <- drop(crossprod(w1, Iy1))
                    }
                    
                    
                    F0_final <- colMeans(F_fold_0)
                    F1_final <- colMeans(F_fold_1)
                    
                    q0 <- quantile_lookup(y_grid, F0_final, taus)$quantile
                    q1 <- quantile_lookup(y_grid, F1_final, taus)$quantile
                    tau_boot_results <- q1 - q0
                    
                    Test_statistics_boot <- sqrt(n)*max(F0_final-F1_final-(F0_hat-F1_hat))
                    
                    
                    
                    return(c(tau_boot_results,Test_statistics_boot))
                    
                  }
end_time <- proc.time()
elapsed_time <- end_time - start_time
cat(sprintf("  -> Total simulation tim: %.2f seconds.\n", elapsed_time[3]))


#p value
pvalues<-mean(result[,18]>Test_statistics)


#CI and CB
qte_boot <- result[,1:17]
bootstrap_mean<-colMeans(qte_boot)
se <- apply(qte_boot, 2, sd, na.rm = TRUE)

alpha<-0.1

ci_lower <- qte_point-qnorm(1 - alpha / 2)*se
ci_upper <- qte_point+qnorm(1 - alpha / 2)*se

qte_boot_center<-apply(qte_boot, 1, function(x) abs(x-QTE_Epy3)/se)
qte_boot_center_max<-apply(qte_boot_center,1,max)

crtical_values <- quantile(qte_boot_center_max,1-alpha)

cb_lower <- qte_point-crtical_values*se
cb_upper <- qte_point+crtical_values*se

# ================================5. Plot====================================
#windows(15,9)
data_plot<-data.frame(cbind(taus,QTE_Epy3,ci_lower,ci_upper,
                            cb_lower,cb_upper))
library('ggplot2')
ggplot(data_plot, aes(x = taus_new)) +
  scale_x_continuous(expand=c(0.015,0))+
  geom_ribbon(aes(ymin = cb_lower, ymax = cb_upper), alpha = 0.2,linetype = 1) +  
  geom_line(aes(y = QTE_Epy3))+
  geom_line(aes(y = -205.4718),color='black',lty=4)+
  labs(x = "Quantiles", y = "Estimated QTE")  +
  theme_test()+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major.x =element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.text = element_text(size=9, color="black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks.x = element_blank())

stopCluster(cl)

