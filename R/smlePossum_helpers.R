#' Observed-data log-likelihood for the sieve maximum likelihood estimator (SMLE)
#'
#' This function returns the value of the observed-data log-likelihood (equation (#) in Lotspeich et al. (2023+))
#' for a given dataset and parameter values `theta` and `p`.
#
#'
#' @param Y Column name with the outcome
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated covariates 
#' @param X_val Column(s) with the validated covariates 
#' @param Z (Optional) Column(s) with additional error-free covariates 
#' @param Bspline Vector of columns containing the B-spline basis functions 
#' @param comp_dat_val Dataset containing rows for validated subjects' data (a matrix)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Vector of columns in \code{data} that pertain to the covariates in the analysis model.
#' @param theta Parameters for the analysis model (a column vector)
#' @param p B-spline coefficients for the approximated covariate error model (a matrix)
#' @return Scalar value of the function

smle_loglik = function(Y = NULL, offset = NULL, X_unval = NULL, X_val = NULL, Z = NULL, Bspline = NULL, comp_dat_val, comp_dat_unval, theta_pred, theta, p) {
  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsX = !is.null(X_unval)
  ## ------ If unvalidated variable was left blank, assume error-free
  # ----------------------------------------- Determine error setting

  if (errorsX) {
    m = nrow(p)
  }

  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  lambda = comp_dat_val[, offset] * exp(as.numeric((cbind(int = 1, comp_dat_val[, theta_pred]) %*% theta)))
  pY_X = dpois(x = comp_dat_val[, Y], lambda = lambda) 
  log_pY_X = log(pY_X)
  log_pY_X[log_pY_X == -Inf] = 0
  return_loglik = sum(log_pY_X)
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  if (errorsX) {
    ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
    pX = p[comp_dat_val[, "k"], ]
    log_pX = log(pX)
    log_pX[log_pX == -Inf] = 0
    return_loglik = return_loglik + sum(comp_dat_val[, Bspline] * log_pX)
    ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  }
  #################################################################################
  # -------------------------------------------------------- For validated subjects

  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  if (!is.null(offset)) {
    lambda = comp_dat_unval[, offset] * 
      exp(as.numeric((cbind(int = 1, comp_dat_unval[, theta_pred]) %*% theta)))
  } else {
    lambda = exp(as.numeric((cbind(int = 1, comp_dat_unval[, theta_pred]) %*% theta)))
  }
  pY_X = dpois(x = comp_dat_unval[, Y], 
               lambda = lambda) 
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  if (errorsX) {
    ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
    pX = p[comp_dat_unval[, "k"], ]
    ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  } else {
    pX = rep(1, nrow(comp_dat_unval[, ]))
  }
  ################################################################################
  ## Calculate sum of P(y|xk) x Bj(X*) x p_kj ------------------------------------
  if (errorsX) {
    person_sum = rowsum(x = pY_X * pX * comp_dat_unval[, Bspline], 
                        group = comp_dat_unval[, "row_num"], 
                        reorder = FALSE)
  }
  person_sum = rowSums(person_sum)
  log_person_sum = log(person_sum)
  log_person_sum[log_person_sum == -Inf] = 0
  ## And sum over them all -------------------------------------------------------
  return_loglik = return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}

#' Profiles out nuisance parameters from the observed-data log-likelihood for a given value of theta
#'
#' For a given vector `theta` to parameterize P(Y|X,Z), this function repeats the EM algorithm to find
#' the values of `gamma` and `p` at convergence. The resulting parameters are used to find the profile
#' log-likelihood for `theta` by plugging them into the observed-data log-likelihood.
#' This function is used by `pl_theta()`.
#
#' @param theta Parameters for the analysis model (a column vector)
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#'
#' @return Profile likelihood for `theta`: the value of the observed-data log-likelihood after profiling out other parameters.
#'
#' @importFrom stats as.formula
#' @importFrom stats glm
#'
#' @noRd

profile_out = function(theta, Y = NULL, offset = NULL, X_unval = NULL, X_val = NULL, Z = NULL, Bspline = NULL,
                        comp_dat_unval, theta_pred, p0, p_val_num, TOL, MAX_ITER) {
  # Save constants --------------------------------------------------
  N = max(comp_dat_unval[, "row_num"]) ## total sample size
  n = min(comp_dat_unval[, "row_num"]) - 1 ## phase II sample size
  sn = ncol(p0) ## number of B-spline sieves
  m = nrow(p0) ## number of distinct values of X
  prev_p = p0 ## starting values for the B-spline coefficients
  
  # For the E-step, save static P(Y|X) for unvalidated --------------
  ### P(Y|X) --------------------------------------------------------
  mu_theta = as.numeric(cbind(int = 1, comp_dat_unval[, theta_pred]) %*% theta)
  lambda = exp(mu_theta)
  if (!is.null(offset)) {
    lambda = comp_dat_unval[, offset] * lambda
  }
  pY_X = dpois(x = comp_dat_unval[, Y], 
               lambda = lambda) 
  ### -------------------------------------------------------- P(Y|X)
  
  # Set parameters for algorithm convergence --------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  
  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(X|X*) -------------------------------------------------------
    ### p_kj ----------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    ### multiply by the B-spline terms
    pX = prev_p[rep(seq(1, m), each = (N - n)), ] * comp_dat_unval[, Bspline]
    ### ---------------------------------------------------------- p_kj
    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------
    ### P(Y|X,Z)p_kjB(X*) -------------------------------------------
    psi_num = c(pY_X) * pX
    ### Update denominator ------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk) ------------------
    psi_denom = rowsum(x = psi_num,
                       group = comp_dat_unval[, "row_num"])
    #### Then sum over the sn splines -------------------------------
    psi_denom = rowSums(psi_denom)
    #### Avoid NaN resulting from dividing by 0 ---------------------
    psi_denom[psi_denom == 0] = 1
    ### And divide them! --------------------------------------------
    psi_t = psi_num / psi_denom
    ### Update the w_kyi for unvalidated subjects -------------------
    ### by summing across the splines/ columns of psi_t -------------
    w_t = rowSums(psi_t)
    ### ----------------------------- Estimate conditional expectations
    # ---------------------------------------------------------- E Step
    ###################################################################
    
    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N ---------
    new_p_num = p_val_num +
      rowsum(x = psi_t,
             group = rep(x = seq(1, m),
                         each = (N - n)),
             reorder = TRUE)
    new_p = t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv = abs(new_p - prev_p) < TOL
    ## -------------------------------------------------- Update {p_kj}
    # ---------------------------------------------------------- M Step
    ###################################################################
    # Check for convergence -------------------------------------------
    CONVERGED = mean(p_conv) == 1
    
    # Update values for next iteration  -------------------------------
    it = it + 1
    prev_p = new_p
    #  ------------------------------- Update values for next iteration
  }
  
  if(it == MAX_ITER & !CONVERGED) {
    CONVERGED_MSG = "MAX_ITER reached"
    new_p = matrix(data = NA,
                    nrow = nrow(p0),
                    ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG = "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}

#' Profile likelihood for theta, the analysis model parameters
#'
#' This function returns the value of the profile log-likelihood for parameters `theta` of the
#' analysis model P(Y|X,Z) after perturbing element `k` of `theta` by some small amount `h_N`.
#
#' @param k A numeric index between 1 and the dimension of theta for the element of theta to be perturbed
#' @param theta Parameters for the analysis model (a column vector) at convergence, resulting from the EM algorithm
#' @param h_N Size of the small perturbation in `theta[k]`, by default chosen to be `h_N =  N ^ ( - 1 / 2)`
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_val Dataset containing rows for validated subjects' data (a matrix)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return Profile likelihood for `theta` after perturbing element `k` by `h_N`.
#' @noRd

pl_theta = function(k, theta, h_N, Y, offset, X_unval, X_val, Z, Bspline, comp_dat_val, comp_dat_unval,
                     theta_pred, p0 = NULL, p_val_num = NULL, TOL, MAX_ITER) {
  pert = theta
  pert[k] = pert[k] + h_N
  pl_params = profile_out(theta = pert,
                           Y = Y,
                           offset = offset,
                           X_unval = X_unval,
                           X_val = X_val,
                           Z = Z,
                           Bspline = Bspline,
                           comp_dat_unval = comp_dat_unval,
                           theta_pred = theta_pred,
                           p0 = p0,
                           p_val_num = p_val_num,
                           TOL = TOL,
                           MAX_ITER = MAX_ITER)
  if(pl_params$converged) {
    od_loglik_pert = smle_loglik(Y = Y,
                                 offset = offset,
                                 X_unval = X_unval,
                                 X_val = X_val,
                                 Z = Z,
                                 Bspline = Bspline,
                                 comp_dat_val = comp_dat_val,
                                 comp_dat_unval = comp_dat_unval,
                                 theta_pred = theta_pred,
                                 theta = pert,
                                 p = pl_params$p_at_conv)
    
  } else { od_loglik_pert = NA }
  return(od_loglik_pert)
}
