rm(list=ls())
gc()

# ==================1. Environment Setup & Libraries=================================

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
source('sim_utils/SLE.R')
source('sim_utils/energy_balance.R')
source('sim_utils/LDML.R')
source('sim_utils/output.R')


# ============================2. Parallel Setup=================================

detectCores(no_cores)
#no_cores<-80     #7709
no_cores<-60     #7707
#no_cores<-8      #mac
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
  
  # Load core code
  source('R/main_estimation.R')
  
  # Load auxiliary code
  source('sim_utils/SLE.R')
  source('sim_utils/energy_balance.R')
  source('sim_utils/LDML.R')
})

# ======3. Data & Configuration Setup==========

taus <- c(0.25, 0.5, 0.75)

datatest <- read.csv('data_scenario1.csv')
datatest<-datatest[,-1]

colnames(datatest)<-c('y','t',paste("X.",1:3,sep=""),'y0','y1')
taus<-c(0.25,0.5,0.75)
true_qte<-quantile(datatest$y1,taus)-quantile(datatest$y0,taus)
true_qte

plot(density(datatest$y0))
lines(density(datatest$y1),col='red')

main_seed<-20260219
#main_seed<-123456

"Z:\\User\\Documents\\Balance_simulations\\202511"
#"/Users/wuyaqian/Documents/202512"
timestamp<-format(Sys.time(), "%m%d-%H")
save_dir <- paste0("Emp_Sim1_06_", timestamp)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE) 
  cat("---Directory created:", save_dir, "---\n")
}  

# ============================4. Run parallel simulations================================
reps <- 100
n_set<-c(1000,2000,4000)
c_set<-c(2)

cat(sprintf("\n--- Starting %d simulation replications on %d cores...\n", reps, no_cores))

all_est_result<-NULL
for(n_level in n_set){
  cat(sprintf("\n--- Running simulations for n = %d ---\n", n_level))
  registerDoRNG(seed = main_seed + n_level) 
  start_time <- proc.time() 
  result <- foreach(k = 1:reps, .combine = 'rbind',
                      .packages = c('quantreg', 'drf', 'quantregForest', 'osqp','ATE.ncb')) %dopar% {
      
      # --- Data Generation ---
      idx_t1_all <- which(datatest$t == 1);idx_t0_all <- which(datatest$t == 0)
      n1_all <- length(idx_t1_all);n0_all <- length(idx_t0_all)
      p1 <- n1_all / (n1_all+n0_all)
      n1 <- round(n_level * p1)
      n0 <- n_level - n1
      index <- c(
        sample(idx_t1_all, n1, replace = FALSE),
        sample(idx_t0_all, n0, replace = FALSE)
      )
                        
      #index<-sample(1:dim(datatest)[1],n_level,replace = FALSE)
      dataset_now<-datatest[index,]
      y <- dataset_now$y
      D <- as.numeric(dataset_now$t)
      X <- dataset_now[,c(3:4)]
     
      
      y_grid<-sort(unique(y))
      
      
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
      fit_energy <- energy_balance(y, D, X, taus, improved=FALSE,normalize = TRUE)
      qte_energy <- fit_energy$QTE 
      
      # 5. Kernel
      fit_kernel <- kernel_balance(y,D,X,taus,nlam=50)
      qte_kernel <- fit_kernel$QTE 
      
      # 7. Entropy
      fit_entropy <- fun_CCDF_balance_new(y,D,X,taus,fun_obj='entropy', 
                                         num_localCCDF=0,num_moment=1,
                                         moments_normlize = TRUE,alpha=0)
      qte_entropy_moments1 <- fit_entropy$final_results$moment_1$QTE 
      
      
      # 8. Initial & DR
      
      method_ccdf <- 'kernel'
      QTE_initial<-NULL;QTE_DR<-NULL;
      QTE_our_moment0<-NULL;QTE_our_moment1<-NULL
      
      for(c_level in c_set){
        bandwidth_c <- c_level
        fit_ccdf_full <- fun_CCDF_balance_old(y, D, X,taus,method_ccdf = method_ccdf,
                                              fun_obj='entropy', num_localCCDF=2, 
                                              num_moment=c(0,1), 
                                              moments_normlize = TRUE,
                                              est_ccdf_normlize=TRUE,
                                              kernel_type=1, 
                                              bandwidth_c=bandwidth_c,alpha=0
        )
        
        qte_initial<-fit_ccdf_full$ccdf_results$QTE 
        qte_DR<-fit_ccdf_full$ccdf_results$QTE_DR 
        
        if(n_level==1000){
          alpha=0.01
        }else{
          alpha=0
        }
        fit_ccdf_balance <- fun_CCDF_balance_new(y, D, X, taus,K=2, method_ccdf=method_ccdf,
                                                 fun_obj='entropy', num_localCCDF=2, 
                                                 num_moment=c(0,1), 
                                                 moments_normlize = TRUE,
                                                 est_ccdf_normlize=TRUE, 
                                                 kernel_type=1, 
                                                 bandwidth_c=bandwidth_c,alpha=alpha)
        
        qte_our_moment0 <- fit_ccdf_balance$final_results$moment_0$est_QTE_cdf_mean 
        qte_our_moment1 <- fit_ccdf_balance$final_results$moment_1$est_QTE_cdf_mean 
        
        QTE_initial<-c(QTE_initial,qte_initial)
        QTE_DR<-c(QTE_DR,qte_DR)
        QTE_our_moment0<-c(QTE_our_moment0,qte_our_moment0)
        QTE_our_moment1<-c(QTE_our_moment1,qte_our_moment1)

      }

      
      est <- c(qte_benchmark, 
               qte_ipw,qte_sle,
               qte_energy, qte_kernel,
               qte_entropy_moments1,
               QTE_initial,QTE_DR,
               QTE_our_moment0, QTE_our_moment1)

      
      #matrix(est,21,3,byrow = TRUE)
      return(est)
    }
    
    end_time <- proc.time()
    elapsed_time <- end_time - start_time
    cat(sprintf("  -> Total simulation time for n=%d: %.2f seconds.\n", n_level, elapsed_time[3]))
    

    # ==========================5. Save Results=================================
    
    names0 <- c('DIQ', 'IPW-PLE','IPW-SLE', 'Energy','Kernel','Epy1')
    names2 <- c('Initial', 'DR', 'Epy2' ,'Epy3')
    
    names3 <- as.vector(sapply(names2, function(p) {
      paste0(p, "_c", 1:length(c_set))
    }))
    est_names<-rep(c(names0,names3),length(taus))
 
    
    R=output(true_qte,result)
    rownames(R)<-est_names
    R
    
    file_name <- paste0(sprintf("emp_sim1_n=%d_", n_level),timestamp,".csv")
    full_path <- file.path(save_dir, file_name)
    write.csv(R,file = full_path, quote = FALSE)
    cat(sprintf("  -> Saving results to: %s\n", file_name))
    
    all_est_result<-cbind(all_est_result,R,Gap = NA)
}

file_name2 <- paste0("all_emp_sim_",timestamp,".csv")
full_path2 <- file.path(save_dir, file_name2)
write.csv(all_est_result,file = full_path2,na = "")
cat(sprintf("--- Saving all results to: %s\n", file_name2))


parallel::stopCluster(cl)
cat("Cluster stopped.\n")



