library(osqp)

get_balanced_weights<-function(tau,D_val,y_train,D_train,X_train,CCDF_train,y_grid,num_localCCDF,num_balance_grid,num_moment,
                               fun_obj = fun_obj,moments_normlize=moments_normlize,delta=delta,...){
  
  n_train <- length(y_train)
  
  design <- build_balance_design(
    tau = tau,
    y_train = y_train,
    D_train = D_train,
    X_train = X_train,
    CCDF_train = CCDF_train,
    y_grid = y_grid,
    num_localCCDF = num_localCCDF,
    num_balance_grid = num_balance_grid,
    num_moment = num_moment,
    moments_normlize = moments_normlize
  )
  
  bal_features_arm <- design$bal_features[D_train == D_val, , drop = FALSE]
  n_arm <- nrow(bal_features_arm)
  if (n_arm == 0) {
    stop("No observations in selected treatment arm for weighting.")
  }
  
  if(fun_obj=='entropy'){
    base_weights <-  rep(1, n_arm)
  }else if(fun_obj=='variance'){
    base_weights <-  rep(n_train/n_arm, n_arm)
  }
  
  
  fit_w<-weight_w(bal_features_arm, design$target_means, base_weights, obj = fun_obj,delta=delta)
  
  
  balanced_weights <- fit_w$w / sum(fit_w$w)
  
  return(balanced_weights=balanced_weights)
}


build_balance_design <- function(tau, y_train, D_train, X_train, CCDF_train, y_grid,
                                 num_localCCDF, num_balance_grid, num_moment,
                                 moments_normlize = FALSE) {
  marginal_CDF_train <- colMeans(CCDF_train, na.rm = TRUE)
  marginal_CDF_train_scale <- marginal_CDF_train / max(marginal_CDF_train)
  marginal_CDF_train <- as.numeric(
    lapply(y_grid, function(x) max(0, max(marginal_CDF_train_scale[y_grid <= x])))
  )
  
  q_pre <- min(y_grid[marginal_CDF_train >= tau])
  
  # 1. Construct local balancing features (based on CCDF)
  local_features <- switch(
    as.character(num_localCCDF),
    "0" = NULL,
    "1" = {
      matrix(CCDF_train[, which(y_grid == q_pre)], ncol = 1)
    },
    "2" = {
      tau_sets <- tau + c(-0.05,-0.03, 0,0.03, 0.05)
      tau_sets <- .clamp_prob(tau_sets)
      q_pre_sets <- sapply(tau_sets, function(x) min(y_grid[marginal_CDF_train >= x]))
      q_idx <- sapply(q_pre_sets, function(t) which(y_grid == t))
      CCDF_train[, q_idx, drop = FALSE]
    },
    "3" = {
      q_idx <- floor(seq(1, length(y_grid), length.out = num_balance_grid))
      CCDF_train[, q_idx, drop = FALSE]
    },
    "4" = {
      q_pre_sets1 <- quantile(y_grid, seq(0.1, 0.9, length.out = num_balance_grid))
      q_idx1 <- sapply(q_pre_sets1, function(t) which(y_grid == t))
      tau_sets <- tau + c(-0.05,-0.03, 0,0.03, 0.05)
      tau_sets <- .clamp_prob(tau_sets)
      q_pre_sets2 <- sapply(tau_sets, function(x) min(y_grid[marginal_CDF_train >= x]))
      q_idx2 <- sapply(q_pre_sets2, function(t) which(y_grid == t))
      q_idx <- sort(unique(c(q_idx1, q_idx2)))
      CCDF_train[, q_idx, drop = FALSE]
    },
    stop("Invalid num_localCCDF value.")
    
  )
  
  # 2. Construct moment balancing features
  X_train_trans <- if (moments_normlize) scale(X_train) else X_train
  
  moment_features <- switch(
    as.character(num_moment),
    "0" = NULL,
    "1" = X_train_trans,
    "2" = cbind(X_train_trans, X_train_trans^2),
    stop("Invalid num_moment value.")
  )
  
  
  # 3. Combine features and calculate weights
  
  bal_features <- cbind(local_features, moment_features)
  
  if (!is.null(bal_features) && ncol(bal_features) > 0) {
    bal_features <- as.matrix(bal_features)
    for (j in seq_len(ncol(bal_features))) {
      col <- bal_features[, j]
      if (!all(is.finite(col))) {
        col_mean <- mean(col[is.finite(col)], na.rm = TRUE)
        if (!is.finite(col_mean)) col_mean <- 0
        col[!is.finite(col)] <- col_mean
        bal_features[, j] <- col
      }
    }
  }
  bal_features <- .drop_zero_var_cols(bal_features)
  
  if (ncol(bal_features) == 0) {
    target_means <- numeric(0)
  } else {
    target_means <- colMeans(bal_features, na.rm = TRUE)
    finite_col <- is.finite(target_means)
    if (!all(finite_col)) {
      bal_features <- bal_features[, finite_col, drop = FALSE]
      target_means <- target_means[finite_col]
    }
  }
  
  return(list(
    bal_features = bal_features,
    target_means = target_means
  ))
}

weight_w<- function(C, M, Q,obj, delta) {
  #C : balance condition
  #M : target
  #Returns weights for control group
  #Q: based weight
  if(missing(obj)){obj='entropy'}
  
  if(obj=='entropy'){
    w_fit<-weight_entropy(C, M, Q, delta)
  }else if (obj=='variance'){
    w_fit<-weight_var(C, M, Q)
  }
  return(list(Z = w_fit$Z,
              w = w_fit$w,
              opt.out = w_fit$opt.out,
              convergence=w_fit$convergence))
}

weight_entropy <-  function(C, M, Q, delta) {
  #C : balance condition
  #M : target
  #Q: based weight
  
  C<-as.matrix(C)
  nc <- nrow(C)
  
  #Z: lagrangian multiplier
  #W: not the final weights,just a  function
  W <- function(Z) {
    drop((1/Q) *exp(-C %*% Z))
  }
  
  objective.EB <- function(Z) {
    log(sum(W(Z))) + sum(M * Z) + delta * sum(Z^2)
  }
  
  gradient.EB <- function(Z) {
    w <- W(Z)
    drop(M - w %*% C/sum(w)+ 2 * delta * Z)
  }
  
  opt.out <- optim(par = rep(0, ncol(C)),
                   fn = objective.EB,
                   gr = gradient.EB,
                   method = "BFGS",
                   control = list(trace = 0,
                                  reltol =  sqrt(.Machine$double.eps),
                                  maxit = 1000))
  
  w <- W(opt.out$par)/sum(W(opt.out$par))
  
  return(list(Z = setNames(opt.out$par, colnames(C)),
              w = w,
              opt.out = opt.out,
              convergence=opt.out$convergence))
}


weight_var <-  function(C, M,Q) {
  #C : balance condition
  #M : target
  #Returns weights for control group
  #Q: based weight
  
  C<-as.matrix(C)
  nc<-dim(C)[1]
  In<-diag(nc)
  
  Amat<-rbind(diag(nc),t(C),1)
  
  min.w = 1e-8;
  lvec<-c(rep(min.w ,nc),M,1)
  uvec<-c(rep(Inf,nc),M,1)
  
  options.list<-osqpSettings(verbose=FALSE)
  
  fit<-solve_osqp(P=In,q=-Q,A=Amat,l=lvec,u=uvec,options.list)
  w<-fit$x
  
  return(list(w = w))
}

.clamp_prob <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

.uniform_weights <- function(n) {
  rep(1 / n, n)
}

.drop_zero_var_cols <- function(M, tol = 1e-12) {
  if (is.null(M)) return(NULL)
  M <- as.matrix(M)
  if (ncol(M) == 0) return(M)
  keep <- vapply(seq_len(ncol(M)), function(j) {
    col <- M[, j]
    finite <- is.finite(col)
    if (sum(finite) <= 1) return(FALSE)
    stats::var(col[finite]) > tol
  }, logical(1))
  if (!any(keep)) {
    return(matrix(numeric(0), nrow = nrow(M), ncol = 0))
  }
  M[, keep, drop = FALSE]
}

.context_str <- function(context = NULL) {
  if (is.null(context)) return("")
  keys <- c("fold", "tau", "D_val", "moment")
  vals <- vapply(keys, function(k) {
    if (!is.null(context[[k]])) as.character(context[[k]]) else NA_character_
  }, character(1))
  keep <- !is.na(vals)
  if (!any(keep)) return("")
  paste0(" [", paste0(keys[keep], "=", vals[keep], collapse = ", "), "]")
}
