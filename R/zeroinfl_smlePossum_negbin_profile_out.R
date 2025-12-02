smlePossum_negbin_profile_out = function(beta, eta, theta, N, n, beta_cols, eta_cols, Y, Bspline, comp_dat_all,
                                         p0, p_val_num, tol, max_iter) {
  # Save useful constants -------------------------------------------
  ## Dimensions and starting values ---------------------------------
  sn = ncol(p0)
  m = nrow(p0)
  prev_p = p0

  ## Create design matrix for P(Y|X,C) model ------------------------
  ### Convert complete data to matrix -------------------------------
  comp_dat_all = as.matrix(comp_dat_all)

  ## Split off complete data for unvalidated rows -------------------
  comp_dat_unval = comp_dat_all[-c(1:n), ]

  # Estimate p using EM -----------------------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  while(it <= max_iter & !CONVERGED) {
    ############################################################################
    # E Step -------------------------------------------------------------------
    E_step_res = E_step_zeroinfl_nb(beta_cols = beta_cols,
                                    eta_cols = eta_cols,
                                    Y = Y,
                                    beta = beta,
                                    eta = eta,
                                    theta = theta,
                                    prev_p = prev_p,
                                    Bspline = Bspline,
                                    comp_dat_unval = comp_dat_unval,
                                    m = m,
                                    N = N,
                                    n = n,
                                    use_predict_pYgivX = FALSE) ### use built-in PMF with values provided
    ############################################################################
    # M Step (but only update the p_{kj}) --------------------------------------
    M_step_res = M_step_zeroinfl_nb_ponly(phi_aug = E_step_res$phi_aug,
                                          psi_t = E_step_res$psi_t,
                                          prev_p = prev_p,
                                          p_val_num = p_val_num,
                                          m = m,
                                          N = N,
                                          n = n,
                                          tol = tol)
    ############################################################################
    # Check for global convergence ---------------------------------------------
    CONVERGED = M_step_res$prop_conv == 1
    # Update values for next iteration  ----------------------------------------
    it = it + 1
    prev_p = M_step_res$new_p
  }

  if(it > max_iter & !CONVERGED) {
    CONVERGED_MSG = "max_iter reached"
    new_p = matrix(data = NA,
                    nrow = nrow(p0),
                    ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG = "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = E_step_res$psi_t,
              "p_at_conv" = prev_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
