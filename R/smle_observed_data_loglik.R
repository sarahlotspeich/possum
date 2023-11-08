#' Observed-data log-likelihood for the sieve maximum likelihood estimator (SMLE)
#'
#' This function returns the value of the observed-data log-likelihood (equation (#) in Lotspeich et al. (2023+))
#' for a given dataset and parameter values `theta` and `p`.
#
#'
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y Column name with the outcome
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offse
#' @param X_unval Column(s) with the unvalidated covariates 
#' @param X_val Column(s) with the validated covariates 
#' @param Z (Optional) Column(s) with additional error-free covariates 
#' @param Bspline Vector of columns containing the B-spline basis functions 
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Vector of columns in \code{data} that pertain to the covariates in the analysis model.
#' @param theta Parameters for the analysis model (a column vector)
#' @param p B-spline coefficients for the approximated covariate error model (a matrix)
#' @return Scalar value of the function

smle_observed_data_loglik = function(N, n, Y = NULL, offset = NULL, X_unval = NULL, X_val = NULL, Z = NULL, Bspline = NULL, comp_dat_all, theta_pred, gamma_pred, theta, gamma, p) {
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
  lambda = comp_dat_all[c(1:n), offset] * exp(as.numeric((cbind(int = 1, comp_dat_all[c(1:n), theta_pred]) %*% theta)))
  pY_X = lambda ^ as.vector(comp_dat_all[c(1:n), Y]) * exp(- lambda) / as.vector(factorial(comp_dat_all[c(1:n), Y]))
  return_loglik = sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  if (errorsX) {
    ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
    pX = p[comp_dat_all[c(1:n), "k"], ]
    log_pX = log(pX)
    log_pX[log_pX == -Inf] = 0
    return_loglik = return_loglik + sum(comp_dat_all[c(1:n), Bspline] * log_pX)
    ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  }
  #################################################################################
  # -------------------------------------------------------- For validated subjects

  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  lambda = comp_dat_all[-c(1:n), offset] * exp(as.numeric((cbind(int = 1, comp_dat_all[-c(1:n), theta_pred]) %*% theta)))
  pY_X = lambda ^ as.vector(comp_dat_all[-c(1:n), Y]) * exp(- lambda) / as.vector(factorial(comp_dat_all[-c(1:n), Y]))
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  if (errorsX) {
    ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
    pX = p[comp_dat_all[-c(1:n), "k"], ]
    ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  } else {
    pX = rep(1, nrow(comp_dat_all[-c(1:n), ]))
  }
  ################################################################################
  ## Calculate sum of P(y|xk) x Bj(X*) x p_kj ------------------------------------
  if (errorsX) {
    person_sum = rowsum(x = c(pY_X * pX) * comp_dat_all[-c(1:n), Bspline], 
                        group = rep(seq(1, (N - n)), times = m), 
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
