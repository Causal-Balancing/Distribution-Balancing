rm(list=ls())
gc()

# ================1. Environment Setup & Libraries===========================

# .libPaths(...) removed for portability


# Load necessary packages
packages.needed <- c('foreach', 'parallel', 'doParallel', 'doRNG', 'quantreg', 'drf', 
                     'quantregForest', 'osqp','ATE.ncb','BMisc')
#new.packages <- packages.needed[!(packages.needed %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
lapply(packages.needed, library, character.only = TRUE)


# setwd(...) removed for portability; run from repository root
# Load core functions (assuming current working dir is the project root)
source('R/main_estimation.R')
# Load auxiliary scripts (Data generation, benchmarks)
source('sim_utils/sim_datagen.R')
source('sim_utils/SLE.R')
source('sim_utils/energy_balance.R')
source('sim_utils/LDML.R')
source('sim_utils/output.R')


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
  library(ATE.ncb)
  source('R/main_estimation.R')
  # Load auxiliary code
  source('sim_utils/sim_datagen.R')
  source('sim_utils/SLE.R')
  source('sim_utils/energy_balance.R')
  source('sim_utils/LDML.R')
})


# =====3. DGP & Create Directory to Store Simulation Results===============

#Calculate True Value (using large sample)

taus <- c(0.25,0.9)

set.seed(main_seed);
current_dgp_func <- data1
dgp_name_str <- "data1"
datatest <- current_dgp_func(1000000)
true_QTE <- quantile(datatest$y1, taus) - quantile(datatest$y0, taus)

true_QTE 

plot(density(datatest$y1),ylim=c(0,1));#lines(density(datatest$y0),col='red')
plot(density(datatest$y[datatest$t==1]),ylim=c(0,1));#lines(density(datatest$y[datatest$t==0]),col='red')

#Create Directory to Store Simulation Results

# output directory is created relative to current working directory
timestamp<-format(Sys.time(), "%m%d-%H")
save_dir <- paste0("B_",dgp_name_str,"_", timestamp)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE) 
  cat("---Directory created:", save_dir, "---\n")
}  

# ================= 4. Run parallel simulations===================================
reps <- 500
n_set<-c(1000,2000,4000)
c_set<-c(1,2,3)

cat(sprintf("\n--- Starting %d simulation replications on %d cores...\n", reps, no_cores))

all_est_result<-NULL
for(n_level in n_set){
  
  cat(sprintf("\n--- Running simulations for n = %d ---\n", n_level))
  registerDoRNG(seed = main_seed + n_level) 
  
  start_time <- proc.time() 
  result <- foreach(k = 1:reps, .combine = 'rbind',
                      .packages = c('quantreg', 'drf', 'quantregForest', 'osqp','ATE.ncb')) %dopar% {
      
      # --- Data Generation ---
      data <- current_dgp_func(n_level)
      X <- data[, 2:3, drop = FALSE]
      y <- data$y
      D <- data$t
      y_grid<-sort(unique(y))
      
      #plot(density(X[D==0,1]));lines(density(X[D==1,1]),col='red')
      #plot(density(X[D==0,2]));lines(density(X[D==1,2]),col='red')
      # --- Estimation ---
      
      # 1.Benchmark
      qte_benchmark <- quantile(y[D==1], taus) - quantile(y[D==0], taus)
      
      # 2. IPW
      fit_ipw <- ipw_fit(y, D, X, taus)
      qte_ipw <- fit_ipw$QTE
      
      # 3. SLE
      fit_sle <- SLE(y, D, X, 2, taus)
      qte_sle <- fit_sle$QTE
      
      # 4. Energy
      time_energy <- system.time({
        fit_energy1 <- energy_balance(y, D, X, taus, improved=FALSE)
      })
      qte_energy1 <- fit_energy1$QTE
      
      
      # 5. Kernel
      fit_kernel <- kernel_balance(y,D,X,taus,nlam=50)
      qte_kernel <- fit_kernel$QTE
      
      
      # 7. Entropy
      fit_entropy <- fun_CCDF_balance_new(y,D,X,taus,fun_obj='entropy', 
                                         num_localCCDF=0,num_moment=1,moments_normlize = TRUE)
      qte_entropy_moments1 <- fit_entropy$final_results$moment_1$QTE
      
      
      # 8. Initial & DR
      method_ccdf <- 'kernel'
      QTE_initial<-NULL;QTE_DR<-NULL;
      QTE_old_moment0<-NULL;QTE_old_moment1<-NULL;
      QTE_our_moment0<-NULL;QTE_our_moment1<-NULL
      
      for(c_level in c_set){
        bandwidth_c <- c_level
        
      
        fit_ccdf_full <- fun_CCDF_balance_old(y, D, X,taus,method_ccdf = 'kernel',
                                              fun_obj='entropy', num_localCCDF=2,
                                              #num_grid=10,
                                              num_moment=c(0,1), 
                                              moments_normlize = TRUE,
                                              est_ccdf_normlize=FALSE, kernel_type=1, 
                                              bandwidth_c=bandwidth_c)
        
        qte_initial<-fit_ccdf_full$ccdf_results$QTE
        qte_DR<-fit_ccdf_full$ccdf_results$QTE_DR
        
  
        
        # 9. Our method
        #bandwidth_c=0.5
        time_epy_run <- system.time({
          fit_ccdf_balance <- fun_CCDF_balance_new(y, D, X, taus,K=2, method_ccdf='kernel',
                                                 fun_obj='entropy', num_localCCDF=2,
                                                 #num_grid=10,
                                                 num_moment=c(0,1), 
                                                 moments_normlize = TRUE,
                                                 est_ccdf_normlize=FALSE, kernel_type=1, 
                                                 bandwidth_c=bandwidth_c)
        })
        
        qte_our_moment0 <- fit_ccdf_balance$final_results$moment_0$est_QTE_cdf_mean
        qte_our_moment1 <- fit_ccdf_balance$final_results$moment_1$est_QTE_cdf_mean
        
        
        
        QTE_initial<-c(QTE_initial,qte_initial)
        QTE_DR<-c(QTE_DR,qte_DR)
        QTE_our_moment0<-c(QTE_our_moment0,qte_our_moment0)
        QTE_our_moment1<-c(QTE_our_moment1,qte_our_moment1)
      }
      

      
      est <- c(qte_benchmark, 
               qte_ipw, qte_sle, 
               qte_energy1,qte_kernel,
               qte_entropy_moments1,
               QTE_initial,QTE_DR,
               QTE_our_moment0, QTE_our_moment1)
      
      #matrix(est,21,3,byrow = TRUE)
      return(c(est,time_energy,time_epy_run))
    }
    
    end_time <- proc.time()
    elapsed_time <- end_time - start_time
    cat(sprintf("  -> Total simulation time for n=%d: %.2f seconds.\n", n_level, elapsed_time[3]))
    
    # ========================5. Save Results====================================
    
    num_est <- (6 + 4 * length(c_set)) * length(taus)
    
    est_results_only <- result[, 1:num_est]
    time_data <- result[, (num_est + 1):ncol(result)]
    
    avg_time<- colMeans(time_data)
    
    cat(sprintf("     [Time Analysis] Average Time for Energy: %.4f seconds\n", avg_time[3]))
    cat(sprintf("     [Time Analysis] Average Time for Epy2 : %.4f seconds\n", avg_time[6]))
    
    names0 <- c('DIQ', 'PLE-IPW', 'SLE-IPW', 'Energy1','Kernel', 'Epy1')
    names2 <- c('Initial', 'DR','Epy2', 'Epy3')

    names3 <- as.vector(sapply(names2, function(p) {
      paste0(p, "_c", 1:length(c_set))
    }))
    est_names<-rep(c(names0,names3),length(taus))
 
    R=output(true_QTE,est_results_only)
    rownames(R)<-est_names
    R
    
    file_name <- paste0(dgp_name_str,sprintf("_n=%d_", n_level),timestamp,".csv")
    full_path <- file.path(save_dir, file_name)
    write.csv(R,file = full_path, quote = FALSE)
    cat(sprintf("  -> Saving results to: %s\n", file_name))
    
    all_est_result<-cbind(all_est_result,R,Gap = NA)
}


file_name2 <- paste0("all_",dgp_name_str,"_",timestamp,".csv")
full_path2 <- file.path(save_dir, file_name2)
write.csv(all_est_result,file = full_path2,na = "")
cat(sprintf("--- Saving all results to: %s\n", file_name2))


parallel::stopCluster(cl)
cat("Cluster stopped.\n")


