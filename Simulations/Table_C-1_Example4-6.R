rm(list=ls())
gc()

# ================1. Environment Setup & Libraries===========================

# .libPaths(...) removed for portability


# Load necessary packages
packages.needed <- c('foreach', 'parallel', 'doParallel', 'doRNG', 'quantreg', 'drf', 
                     'quantregForest', 'osqp','BMisc')
#new.packages <- packages.needed[!(packages.needed %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
lapply(packages.needed, library, character.only = TRUE)


# setwd(...) removed for portability; run from repository root
# Load core functions (assuming current working dir is the project root)
source('R/main_estimation.R')
source('R/inference_bootstrap.R')
source('R/test_bootstrap.R')
# Load auxiliary scripts (Data generation, benchmarks)
source('sim_utils/sim_datagen.R')
source('sim_utils/SLE.R')
source('sim_utils/Other_SD_test.R')



# ======================2. Parallel Setup===================================

available_cores <- parallel::detectCores()
no_cores<-22     #7706
no_cores<-90     #7709
no_cores<-60     #7707

main_seed <- 20260215
cl = parallel::makeCluster(no_cores, type = "SOCK")
registerDoParallel(cl)

# Set random seed stream for reproducibility
parallel::clusterSetRNGStream(cl,iseed = main_seed)

# Export Environment to Cluster
cat("--- Exporting scripts and libraries to cluster workers... ---\n")

# Pre-load code on all workers to avoid repeated sourcing inside loop
clusterEvalQ(cl, {
  library(quantreg)
  library(BMisc)
  library(quantregForest)
  library(drf)
  library(osqp)
  source('R/main_estimation.R')
  source('R/inference_bootstrap.R')
  source('R/test_bootstrap.R')
  # Load auxiliary code
  source('sim_utils/sim_datagen.R')
  source('sim_utils/SLE.R')
})


# =====3. DGP & Create Directory to Store Simulation Results===============

#Calculate True Value (using large sample)

taus <- c(0.5)


# ================= 4. Run parallel simulations===================================
reps <- 1000
num_boot<-999
n_set<-c(1000,2000,4000)
dgp_name_str <- c("data1","data_ks_1","data_ks_1")
alpha<-c(0.05,0.1)

cat(sprintf("\n--- Starting %d simulation replications on %d cores...\n", reps, no_cores))

for (dgp_name in dgp_name_str) {
  cat(sprintf("\n--- Running simulations for DGP = %s ---\n", dgp_name))
  for(n_level in n_set){
    
    cat(sprintf("\n--- Running simulations for n = %d ---\n", n_level))
    registerDoRNG(seed = main_seed + n_level) 
    
    start_time <- proc.time() 
    result <- foreach(k = 1:reps, .combine = 'rbind',
                      .packages = c('quantreg', 'drf', 'quantregForest', 'osqp','ATE.ncb')) %dopar% {
                        
                        # --- Data Generation ---
                        data <- match.fun(dgp_name)(n_level)
                        X <- data$X
                        y <- data$y
                        D <- data$t
                        y_grid<-sort(unique(y))
                        
                        if(dgp_name=='data1'){bandwidth_c=1
                        }else{
                          bandwidth_c=0.5
                        }
                
                        
                        fit_inference<-fun_CCDF_balance_with_test(y, D, X, taus, K = 2, fun_obj='entropy', 
                                                                  method_ccdf='kernel',kernel_type=1, bandwidth_c=bandwidth_c,
                                                                  num_localCCDF=3,est_ccdf_normlize=FALSE, num_balance_grid = 50,
                                                                  num_moment=c(0,1), moments_normlize = TRUE,
                                                                  num_boot = num_boot, alpha = 0.05,delta = 0.01,
                                                                  bootstrap_dist = "exponential", seed = NULL,
                                                                  maxit = 500 )
                        
                        pvalue_Epy2<-fit_inference$test$moment_0$pvalue
                        pvalue_Epy3<-fit_inference$test$moment_1$pvalue
                        
                        
                        rejection_Epy2 <- c(as.numeric(pvalue_Epy2<alpha[1]),as.numeric(pvalue_Epy2<alpha[2]))
                        rejection_Epy3 <- c(as.numeric(pvalue_Epy3<alpha[1]),as.numeric(pvalue_Epy3<alpha[2]))
                        
                        
                        fit_DH_test<-DH_test(y,D,X,d=2,B=num_boot,alpha=alpha)
                        
                        rejection_DH<-fit_DH_test$rejection
                        
                        fit_SLE_test<-SLE_test(y,D,X,d=2,B=num_boot,alpha=alpha)
                        
                        rejection_SLE<-fit_SLE_test$rejection
                        
                        est <- c(rejection_SLE,rejection_DH,rejection_Epy2,rejection_Epy3)
                        return(est)
                      }
    
    end_time <- proc.time()
    elapsed_time <- end_time - start_time
    cat(sprintf("  -> Total simulation time for n=%d: %.2f seconds.\n", n_level, elapsed_time[3]))
    
    # ========================5. Save Results====================================
    
    est_names <- rep(c('rejection_SLE','rejection_DH','rejection_Epy2','rejection_Epy3'), 
                     each = length(alpha))
    
    colnames(result) <- est_names
    
    
    avg_rej_rate <- colMeans(result, na.rm = TRUE)
    
    
    cat("     [Rejection Rate] 'SLE','SLE','DH','DH','Epy2','Epy2','Epy3','Epy3':\n")
    print(avg_rej_rate)
    
  }
  
}



stopCluster(cl)



