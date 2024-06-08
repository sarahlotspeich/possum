#' Observed-data log-likelihood for the sieve maximum likelihood estimator (SMLE)
#'
#' This function returns the value of the observed-data log-likelihood (equation (#) in Lotspeich et al. (2023+))
#' for a given dataset and parameter values `beta` and `p`.
#
#'
#' @param Y Column name with the outcome
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated covariates 
#' @param X Column(s) with the validated covariates 
#' @param Z (Optional) Column(s) with additional error-free covariates 
#' @param comp_dat_val Dataset containing rows for validated subjects' data (a matrix)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param beta_pred Vector of columns in \code{data} that pertain to the covariates in the analysis model.
#' @param beta Parameters for the analysis model (a column vector)
#' @param ppv Positive predictive value for the misclassification mechanism (a scalar)
#' @return Scalar value of the function

smle_loglik = function(Y = NULL, offset = NULL, X_unval = NULL, X = NULL, Z = NULL, comp_dat_val, comp_dat_unval, beta_pred, beta, ppv) {
  ##############################################################################
  # Save useful constants ------------------------------------------------------
  N = max(comp_dat_unval[, "row_num"]) ## total sample size 
  n = max(comp_dat_val[, "row_num"]) ## validation sub-sample size
  
  # Create combined dataset of validated and unvalidated -----------------------
  comp_dat_all = rbind(comp_dat_val, 
                       comp_dat_unval)
  ##############################################################################
  # Calculate probabilities ----------------------------------------------------
  ## Analysis model: P(Y|X) ----------------------------------------------------
  ### mu = beta0 + beta1X + beta2Z + ... 
  mu_beta = as.numeric(cbind(int = 1, comp_dat_all[, beta_pred]) %*% beta)
  ### lambda = exp(beta0 + beta1X + beta2Z + ... )
  lambda = exp(mu_beta)
  ### If offset specified, lambda = offset x exp(beta0 + beta1X + beta2Z + ... )
  if (!is.null(offset)) {
    lambda = comp_dat_all[, offset] * lambda
  }
  ### Calculate P(Y|X) from Poisson distribution 
  pYgivX = dpois(x = comp_dat_all[, Y], 
                 lambda = lambda)
  ##############################################################################
  ## Misclassification mechanism: P(X|X*) --------------------------------------
  ### P(X|X*=1) = ppv ^ X * (1 - X) 
  pXgivXstar = ppv ^ comp_dat_all[, X] * (1 - ppv) ^ (1 - comp_dat_all[, X])
  #### P(X|X*=0) = I(X=0) since false negatives aren't possible 
  pXgivXstar[comp_dat_all[, X_unval] == 0] = as.numeric(comp_dat_all[comp_dat_all[, X_unval] == 0, X] == 0)
  ##############################################################################
  ## Joint conditional: P(Y,X|X*) ----------------------------------------------
  pYXgivXstar = pYgivX * pXgivXstar
  ##############################################################################
  # Log-likelihood contribution of validated observations  ---------------------
  ## Sum over log{P(Y|X)} for validated observations (first n rows)
  log_pYgivX = log(pYgivX[1:n])
  log_pYgivX[log_pYgivX == -Inf] = 0
  return_loglik = sum(log_pYgivX)
  ## Sum over log{P(X|X*)} 
  log_pXgivXstar = log(pXgivXstar[1:n])
  log_pXgivXstar[log_pXgivXstar == -Inf] = 0
  return_loglik = return_loglik + sum(log_pXgivXstar)
  ##############################################################################
  # Log-likelihood contribution of unvalidated observations  -------------------
  ## Marginalize out X: P(Y|X*) ------------------------------------------------
  pYXgivXstar_unval = pYXgivXstar[-c(1:n)] ### remove validated observations (first n rows)
  pYgivXstar_unval = pYXgivXstar_unval[1:(N - n)] + ### sum over P(Y,X=0|X*) (first N - n rows)
    pYXgivXstar_unval[-c(1:(N - n))] ### and P(Y,X=1|X*) (last N - n rows)
  ## Sum over log{P(Y|X*)}
  log_pYgivXstar_unval = log(pYgivXstar_unval)
  log_pYgivXstar_unval[log_pYgivXstar_unval == -Inf] = 0
  return_loglik = return_loglik + sum(log_pYgivXstar_unval)
  return(return_loglik)
}

#' Profile likelihood for beta, the analysis model parameters
#'
#' This function returns the value of the profile log-likelihood for parameters `beta` of the
#' analysis model P(Y|X,Z) after perturbing element `k` of `beta` by some small amount `h_N`.
#
#' @param k A numeric index between 1 and the dimension of beta for the element of beta to be perturbed
#' @param beta Parameters for the analysis model (a column vector) at convergence, resulting from the EM algorithm
#' @param h_N Size of the small perturbation in `beta[k]`, by default chosen to be `h_N =  N ^ ( - 1 / 2)`
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param comp_dat_val Dataset containing rows for validated subjects' data (a matrix)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param beta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return Profile likelihood for `beta` after perturbing element `k` by `h_N`.
#' @noRd

pl_beta = function(k, beta, h_N, Y, offset, X_unval, X, Z, comp_dat_val, comp_dat_unval, beta_pred, ppv0, ppv_num0, ppv_denom0, TOL, MAX_ITER) {
  # Define perturbed beta vector -----------------------------------------------
  pert = beta
  pert[k] = pert[k] + h_N
  
  # Save constants -------------------------------------------------------------
  n = max(comp_dat_val[, "row_num"]) ## validation sub-sample size
  N = max(comp_dat_unval[, "row_num"]) ## total sample size
  
  # Estimate *just* ppv using EM algorithm -------------------------------------
  ## Set parameters for algorithm convergence ----------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  ## Initialize PPV ------------------------------------------------------------
  prev_ppv = ppv0
  ## Analysis model: P(Y|X) ----------------------------------------------------
  ### mu = beta0 + beta1X + beta2Z + ... 
  mu_beta = as.numeric(cbind(int = 1, comp_dat_unval[, beta_pred]) %*% beta)
  ### lambda = exp(beta0 + beta1X + beta2Z + ... )
  lambda = exp(mu_beta)
  ### If offset specified, lambda = offset x exp(beta0 + beta1X + beta2Z + ... )
  if (!is.null(offset)) {
    lambda = comp_dat_unval[, offset] * lambda
  }
  ### Calculate P(Y|X) from Poisson distribution -------------------------------
  pYgivX = dpois(x = comp_dat_unval[, Y], 
                 lambda = lambda)
  ##############################################################################
  ## Begin algorithm -----------------------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step -------------------------------------------------------------------
    ## Update the psi_xi = P(X=x|Yi,Xi*,Z) for unvalidated subjects ------------
    ### Misclassification mechanism: P(X|X*) -----------------------------------
    #### P(X|X*=1) = ppv ^ X * (1 - X) -----------------------------------------
    pXgivXstar = prev_ppv ^ comp_dat_unval[, X] * (1 - prev_ppv) ^ (1 - comp_dat_unval[, X])
    #### P(X|X*=0) = I(X=0) since false negatives aren't possible --------------
    pXgivXstar[comp_dat_unval[, X_unval] == 0] = as.numeric(comp_dat_unval[comp_dat_unval[, X_unval] == 0, X] == 0)
    ############################################################################
    ## Estimate conditional expectations --------------------------------------
    ### Update numerator -------------------------------------------------------
    #### P(Y|X,Z)P(X|X*) -------------------------------------------------------
    psi_num = pYgivX * pXgivXstar ##### dim: 2(N - n) x 1
    psi_num_wide = matrix(data = psi_num, 
                          nrow = (N - n), 
                          ncol = 2, 
                          byrow = FALSE)
    ### Update denominator -----------------------------------------------------
    #### P(Y|X=0,Z)P(X=0|X*) + P(Y|X=1,Z)P(X=1|X*) -----------------------------
    psi_denom = rowSums(psi_num_wide) ##### dim: (N - n) x 1
    #### Avoid NaN resulting from dividing by 0 --------------------------------
    psi_denom[psi_denom == 0] = 1
    ### Divide them to get psi = E{I(X=x)|Y,X*} --------------------------------
    psi = psi_num / rep(x = psi_denom, times = 2) 
    #### Add indicators for validated rows -------------------------------------
    psi_aug = c(rep(x = 1, times = n), psi)
    ############################################################################
    # M Step -------------------------------------------------------------------
    ## Update ppv --------------------------------------------------------------
    ### Update numerator by adding sum over X* x psi_1i for unvalidated --------
    ppv_num = ppv_num0 + 
      sum(comp_dat_unval[-c(1:(N-n)), X_unval] * psi[-c(1:(N-n))])
    ### Divide them ------------------------------------------------------------
    new_ppv = ppv_num / ppv_denom0
    ### Check for ppv convergence ----------------------------------------------
    CONVERGED = abs(new_ppv - prev_ppv) < TOL
    ############################################################################
    # Update values for next iteration  ----------------------------------------
    it = it + 1
    prev_ppv = new_ppv
  }
  if(CONVERGED) {
    ### Save PPV at convergence -------------------------------------------------
    ppv_at_conv = new_ppv
    od_loglik_pert = smle_loglik(Y = Y,
                                 offset = offset,                                          
                                 X_unval = X_unval,
                                 X = X,
                                 Z = Z,
                                 comp_dat_val = comp_dat_val,
                                 comp_dat_unval = comp_dat_unval,
                                 beta_pred = beta_pred,
                                 beta = pert,
                                 ppv = ppv_at_conv)
  } else { od_loglik_pert = NA }
  return(od_loglik_pert)
}
