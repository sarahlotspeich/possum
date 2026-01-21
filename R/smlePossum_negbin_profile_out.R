smlePossum_negbin_profile_out = function(beta, theta, N, n, Y, beta_cols, Bspline, comp_dat_all,
                       p0, p_val_num, tol, max_iter) {
  ##############################################################################
  # Save useful constants ------------------------------------------------------
  ## Dimensions and starting values --------------------------------------------
  sn = ncol(p0)
  m = nrow(p0)
  prev_p = p0
  ## Make sure the beta coefficients are a column vector -----------------------
  beta = matrix(data = beta,
                ncol = 1)
  ## Create design matrix for P(Y|X,C) model -----------------------------------
  ### Only among unvalidated rows ----------------------------------------------
  theta_design_mat = comp_dat_all[-c(1:n), c(beta_cols)]
  ### Convert complete data to matrix ------------------------------------------
  comp_dat_all = as.matrix(comp_dat_all)
  ## Split off complete data for unvalidated rows ------------------------------
  comp_dat_unval = comp_dat_all[-c(1:n), ]
  ## Calculate P(Y|X) for theta, since it won't update -------------------------
  ### Only among unvalidated rows ----------------------------------------------
  #### mu = exp(beta0 + beta1X + beta2Z + ) ... --------------------------------
  mu_beta = exp(as.numeric(theta_design_mat %*% beta))
  #### Calculate P(Y|X,Z) from negative binomial distribution ------------------
  pYgivX = dnbinom(x = comp_dat_all[-c(1:n), Y],
                   size = theta,
                   prob = (theta / (mu_beta + theta)))
  ##############################################################################
  # Estimate p using EM --------------------------------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  while(it <= max_iter & !CONVERGED) {
    # E Step -------------------------------------------------------------------
    E_step_res = E_step_nb(prev_beta = beta,
                           prev_theta = theta,
                           Y = Y,
                           beta_cols = beta_cols,
                           prev_p = prev_p,
                           Bspline = Bspline,
                           comp_dat_unval = comp_dat_unval,
                           m = m,
                           N = N,
                           n = n)
    ############################################################################
    # M Step (but only update the p_{kj}) --------------------------------------
    M_step_res = M_step_nb_ponly(psi_t = E_step_res$psi_t,
                                 prev_p = prev_p,
                                 p_val_num = p_val_num,
                                 m = m,
                                 N = N,
                                 n = n,
                                 tol = tol)
    ############################################################################
    # Check for convergence ----------------------------------------------------
    CONVERGED = M_step_res$prop_conv == 1
    # Update values for next iteration  ----------------------------------------
    it = it + 1
    prev_p = M_step_res$new_p
  }
  # -------------------------------------------------------- Estimate p using EM
  ##############################################################################
  # Check for convergence ------------------------------------------------------  
  if(it > max_iter & !CONVERGED) {
    CONVERGED_MSG = "max_iter reached"
    new_p = matrix(data = NA,
                    nrow = nrow(p0),
                    ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG = "converged"
  # ------------------------------------------------------ Check for convergence 
  return(list("psi_at_conv" = E_step_res$psi_t,
              "p_at_conv" = prev_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
