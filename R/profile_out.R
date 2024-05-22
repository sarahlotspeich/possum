#' Profiles out nuisance parameters from the observed-data log-likelihood for a given value of theta
#'
#' For a given vector `theta` to parameterize P(Y|X,Z), this function repeats the EM algorithm to find
#' the values of `gamma` and `p` at convergence. The resulting parameters are used to find the profile
#' log-likelihood for `theta` by plugging them into the observed-data log-likelihood.
#' This function is used by `pl_theta()`.
#
#' @param theta Parameters for the analysis model (a column vector)
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
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

profile_out <- function(theta, n, N, Y = NULL, offset = NULL, X_unval = NULL, X_val = NULL, Z = NULL, Bspline = NULL,
                           comp_dat_all, theta_pred, p0, p_val_num, TOL, MAX_ITER) {
  sn <- ncol(p0)
  m <- nrow(p0)
  prev_p <- p0
  
  # Convert to matrices
  theta_design_mat <- as.matrix(cbind(int = 1,
                                      comp_dat_all[-c(1:n), theta_pred]))
  comp_dat_all <- as.matrix(comp_dat_all)

  # Split complete data for unvalidated data
  comp_dat_unval <- comp_dat_all[-c(1:n), ]
  
  # For the E-step, save static P(Y|X) for unvalidated --------------
  ### P(Y|X) --------------------------------------------------------
  mu_theta = as.numeric((theta_design_mat %*% theta))
  lambda = exp(mu_theta)
  if (!is.null(offset)) {
    lambda = comp_dat_all[-c(1:n), offset] * lambda
  }
  pY_X = dpois(x = comp_dat_all[-c(1:n), Y], lambda = lambda) 
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
                       group = rep(seq(1, (N - n)), times = m))
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
    p_conv <- abs(new_p - prev_p) < TOL
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
    CONVERGED_MSG <- "MAX_ITER reached"
    new_p <- matrix(data = NA,
                    nrow = nrow(p0),
                    ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG <- "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
