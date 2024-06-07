
#' Profiles out nuisance parameters from the observed-data log-likelihood for a given value of beta
#'
#' For a given vector `beta` to parameterize P(Y|X,Z), this function repeats the EM algorithm to find
#' the values of `ppv` at convergence. The resulting parameters are used to find the profile
#' log-likelihood for `beta` by plugging them into the observed-data log-likelihood.
#' This function is used by `pl_beta()`.
#
#' @param beta Parameters for the analysis model (a column vector)
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param beta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param ppv_num0 Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#'
#' @return Profile likelihood for `beta`: the value of the observed-data log-likelihood after profiling out other parameters.
#'
#' @importFrom stats as.formula
#' @importFrom stats glm
#'
#' @noRd

profile_out = function(beta, Y = NULL, offset = NULL, X_unval = NULL, X_val = NULL, Z = NULL, 
                       comp_dat_unval, beta_pred, ppv_num0, TOL, MAX_ITER) {
  # Save constants --------------------------------------------------
  N = max(comp_dat_unval[, "row_num"]) ## total sample size
  n = min(comp_dat_unval[, "row_num"]) - 1 ## phase II sample size
  sn = ncol(p0) ## number of B-spline sieves
  m = nrow(p0) ## number of distinct values of X
  prev_p = p0 ## starting values for the B-spline coefficients
  
  # For the E-step, save static P(Y|X) for unvalidated --------------
  ### P(Y|X) --------------------------------------------------------
  mu_beta = as.numeric(cbind(int = 1, comp_dat_unval[, beta_pred]) %*% beta)
  lambda = exp(mu_beta)
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
  
  # Estimate beta using EM -------------------------------------------
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
  # ---------------------------------------------- Estimate beta using EM
  return(list("psi_at_conv" = psi_t,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
