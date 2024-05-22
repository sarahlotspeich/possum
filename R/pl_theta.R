#' Profile likelihood for theta, the analysis model parameters
#'
#' This function returns the value of the profile log-likelihood for parameters `theta` of the
#' analysis model P(Y|X,Z) after perturbing element `k` of `theta` by some small amount `h_N`.
#
#' @param k A numeric index between 1 and the dimension of theta for the element of theta to be perturbed
#' @param theta Parameters for the analysis model (a column vector) at convergence, resulting from the EM algorithm
#' @param h_N Size of the small perturbation in `theta[k]`, by default chosen to be `h_N =  N ^ ( - 1 / 2)`
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
#' @return Profile likelihood for `theta` after perturbing element `k` by `h_N`.
#' @noRd

pl_theta <- function(k, theta, h_N, n, N, Y, offset, X_unval, X_val, Z, Bspline, comp_dat_all,
                        theta_pred, p0 = NULL, p_val_num = NULL, TOL, MAX_ITER) {
  pert <- theta
  pert[k] <- pert[k] + h_N
  pl_params <- profile_out(theta = pert,
                           n = n,
                           N = N,
                           Y = Y,
                           offset = offset,
                           X_unval = X_unval,
                           X_val = X_val,
                           Z = Z,
                           Bspline = Bspline,
                           comp_dat_all = comp_dat_all,
                           theta_pred = theta_pred,
                           p0 = p0,
                           p_val_num = p_val_num,
                           TOL = TOL,
                           MAX_ITER = MAX_ITER)
  if(pl_params$converged) {
    od_loglik_pert <- smle_observed_data_loglik(N = N,
                                           n = n,
                                           Y = Y,
                                           offset = offset,
                                           X_unval = X_unval,
                                           X_val = X_val,
                                           Z = Z,
                                           Bspline = Bspline,
                                           comp_dat_all = comp_dat_all,
                                           theta_pred = theta_pred,
                                           theta = pert,
                                           p = pl_params$p_at_conv)

  } else { od_loglik_pert <- NA }
  return(od_loglik_pert)
}
