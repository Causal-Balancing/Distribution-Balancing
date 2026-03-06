rm(list=ls())
gc()

# ================1. Environment Setup & Libraries===========================

.libPaths("Z:\\User\\Documents\\R\\win-library\\4.5")


# Load necessary packages
packages.needed <- c('foreach', 'parallel', 'doParallel', 'doRNG', 'quantreg', 'drf', 
                     'quantregForest', 'osqp','ATE.ncb','BMisc')
#new.packages <- packages.needed[!(packages.needed %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
lapply(packages.needed, library, character.only = TRUE)


setwd("Z:\\User\\Documents\\Balance_simulations\\202511")
# Load core functions (assuming current working dir is the project root)
source('R/main_estimation.R')
# Load auxiliary scripts (Data generation, benchmarks)
source('sim_utils/datagen_new.R')
source('sim_utils/SLE.R')
source('sim_utils/energy_balance.R')
source('sim_utils/LDML.R')
source('sim_utils/output.R')



# ======================2. Parallel Setup===================================

detectCores(no_cores)
no_cores<-22     #7706
no_cores<-90     #7709
no_cores<-40    #7707

cl = parallel::makeCluster(no_cores, type = "SOCK")
registerDoParallel(cl)

# Set random seed stream for reproducibility
parallel::clusterSetRNGStream(cl,)

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
main_seed <- 20260217
taus <- c(0.25,0.5, 0.75,0.9)

set.seed(main_seed);

dgp_name_str <- "Kallus_DGP"
datatest <- makedata.sim(1000000,d=20,ovlp=3)
true_QTE <- quantile(datatest$y1, taus) - quantile(datatest$y0, taus)

true_QTE 

plot(density(datatest$y1),ylim=c(0,1));lines(density(datatest$y0),col='red')
#plot(ecdf(datatest$y1));lines(ecdf(datatest$y0),col='red')
#plot(density(datatest$y[datatest$t==1]),ylim=c(0,1));#lines(density(datatest$y[datatest$t==0]),col='red')

#Create Directory to Store Simulation Results

"Z:\\User\\Documents\\Balance_simulations\\202511"
timestamp<-format(Sys.time(), "%m%d-%H")
save_dir <- paste0("B_",dgp_name_str,"_500_", timestamp)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE) 
  cat("---Directory created:", save_dir, "---\n")
}  

# ================= 4. Run parallel simulations===================================
reps <- 200
n_set<-c(1000,2000,4000)

cat(sprintf("\n--- Starting %d simulation replications on %d cores...\n", reps, no_cores))

all_est_result<-NULL
for(n_level in n_set){
  
  cat(sprintf("\n--- Running simulations for n = %d ---\n", n_level))
  registerDoRNG(seed = main_seed + n_level) 

  start_time <- proc.time() 
  result <- foreach(k = 1:reps, .combine = 'rbind',
                      .packages = c('quantreg', 'drf', 'quantregForest', 'osqp','ATE.ncb')) %dopar% {
      
      # --- Data Generation ---
      #data <- makedata.sim(n_level,d=20)
      #X <- data[, 1:20, drop = FALSE]
      data <- makedata.sim(n_level,d=20,ovlp=3)
      X <- data[,1:20, drop = FALSE]
      y <- data$y
      D <- as.numeric(data$t)
      
      y_grid<-sort(unique(y))
      
      
      # --- Estimation ---
      
      # 1.Benchmark
      #qte_benchmark <- quantile(y[D==1], taus) - quantile(y[D==0], taus)
      
      # 2. IPW
      #fit_ipw <- ipw_fit(y, D, X, taus)
      #qte_ipw <- fit_ipw$QTE
      
      # 3. SLE
      #fit_sle <- SLE(y, D, X, 2, taus)
      #qte_sle <- fit_sle$QTE
      
      # 4. Energy
      #fit_energy <- energy_balance(y, D, X, taus, improved=FALSE)
      #qte_energy <- fit_energy$QTE
      #fit_energy2 <- energy_balance(y, D, X, taus, improved=FALSE,normalize = TRUE)
      #qte_energy2 <- fit_energy2$QTE
      
      # 5. Kernel
      #fit_kernel <- kernel_balance(y,D,X,taus,nlam=50)
      #qte_kernel <- fit_kernel$QTE
      
      # 6. LDML
      form_x = paste(paste("X.", 1:dim(X)[2], sep=""),collapse="+")
      form_y = "y"
      form_t = "t"
      #forest_option defined in 'LDML.R'
      my_forest_option = list(nodesize=1, ntree=1000, na.action=na.omit, replace=TRUE)
      fit_ldml <- est.quantile.ldml(taus, data, form_x, form_t, form_y, 
                                    method_ipw=forest, option_ipw=my_forest_option, method_prop=forest, option_prop=my_forest_option, method_cdf=forest, option_cdf=forest_option, 
                                    K=5, K_ipw=NULL, 
                                    semiadaptive=FALSE,  
                                    trim=c(0.01,0.99), trim.type='none', normalize=T, q.oracle=NULL, avg.eqn=T,
                                    oracle.density=NULL)
      qte_ldml<-fit_ldml$qte
      
      # 7. Entropy
      fit_entropy <- fun_CCDF_balance_new(y,D,X,taus,fun_obj='entropy', 
                                         num_localCCDF=0,num_moment=1,alpha=0,delta = 0)
      qte_entropy_moments1 <- fit_entropy$final_results$moment_1$QTE
      
      
      # 8. Initial & DR
      method_ccdf <-'quantregforest';
      #method_ccdf='distforest';
        
      fit_ccdf_full <- fun_CCDF_balance_old(y, D, X,taus,method_ccdf = method_ccdf,
                                              fun_obj='entropy', num_localCCDF=2, 
                                              num_moment=c(0,1), 
                                              moments_normlize = FALSE,
                                              est_ccdf_normlize=FALSE,alpha=0,delta=0)
        
      qte_initial<-fit_ccdf_full$ccdf_results$QTE
      qte_DR<-fit_ccdf_full$ccdf_results$QTE_DR

        

      fit_ccdf_balance <- fun_CCDF_balance_new(y, D, X, taus,K=2, method_ccdf=method_ccdf,
                                                 fun_obj='entropy', num_localCCDF=2, 
                                                 num_moment=c(0,1), 
                                                 moments_normlize = FALSE,
                                                 est_ccdf_normlize=FALSE,alpha=0.01,delta = 0.0)
        
      qte_our_moment0 <- fit_ccdf_balance$final_results$moment_0$est_QTE_cdf_mean
      qte_our_moment1 <- fit_ccdf_balance$final_results$moment_1$est_QTE_cdf_mean
        
        
      
      est <- c(#qte_benchmark, 
               #qte_ipw,
               #qte_entropy_moments1,
               #qte_energy, qte_kernel,
               qte_ldml,
               qte_initial,qte_DR,
               qte_our_moment0, qte_our_moment1)
      
      #matrix(est,21,3,byrow = TRUE)
      return(est)
    }
    
    end_time <- proc.time()
    elapsed_time <- end_time - start_time
    cat(sprintf("  -> Total simulation time for n=%d: %.2f seconds.\n", n_level, elapsed_time[3]))
    

    # ============================5. Save Results=============================
    
    names0 <- c('LDML')
    names2 <- c('Initial', 'DR', 'Epy2', 'Epy3')
    est_names<-rep(c(names0,names2),length(taus))
    R=output(true_QTE,result)
    rownames(R)<-est_names
    R
    
    file_name <- paste0(dgp_name_str,sprintf("_n=%d_", n_level),timestamp,".csv")
    full_path <- file.path(save_dir, file_name)
    write.csv(R,file = full_path, quote = FALSE)
    cat(sprintf("  -> Saving results to: %s\n", file_name))
    
    all_est_result<-cbind(all_est_result,R,Gap = NA)
}

file_name2 <- paste0("all_data3_",timestamp,".csv")
full_path2 <- file.path(save_dir, file_name2)
write.csv(all_est_result,file = full_path2,na = "")
cat(sprintf("--- Saving all results to: %s\n", file_name2))


parallel::stopCluster(cl)
cat("Cluster stopped.\n")



