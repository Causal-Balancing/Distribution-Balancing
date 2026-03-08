
output <- function(true_qte, est_results) {

  num_total_reps <- nrow(est_results)
  est_results_complete <- na.omit(est_results)
  num_valid_reps <- nrow(est_results_complete)
  
  if (num_valid_reps < num_total_reps) {
    cat(sprintf("Warning: Removed %d rows with NA values. Calculations are based on %d valid replications.\n", 
                num_total_reps - num_valid_reps, num_valid_reps))
  }
  
 
  if (num_valid_reps == 0) {
    warning("No valid replications found after removing NAs. Returning an empty result.")
    num_estimators <- ncol(est_results) / length(true_qte)
    return(matrix(NA, nrow = num_estimators * length(true_qte), ncol = 2, 
                  dimnames = list(NULL, c('MADE', 'RMSE'))))
  }
  
 
  num_estimators <- ncol(est_results_complete) / length(true_qte)
  true_qte_vec <- rep(true_qte,num_estimators)
  M_true_qte <- matrix(true_qte_vec, nrow = num_valid_reps, ncol = length(true_qte_vec), byrow = TRUE)
  
  mdad_vec <- apply(abs(est_results_complete - M_true_qte), 2, median)
  mdad <- matrix(mdad_vec, nrow = num_estimators, ncol = length(true_qte), byrow = TRUE)
  
  rmse_vec <- sqrt(colMeans((est_results_complete - M_true_qte)^2))
  rmse <- matrix(rmse_vec, nrow = num_estimators, ncol = length(true_qte), byrow = TRUE)
  
  results_by_tau <- list()
  for (i in 1:length(true_qte)) {
    results_by_tau[[i]] <- cbind(mdad[, i], rmse[, i])
  }
  
  R <- do.call(rbind, results_by_tau)
  R <- round(R,4)
  
  colnames(R) <- c('MADE', 'RMSE')
  
  return(R)
}
