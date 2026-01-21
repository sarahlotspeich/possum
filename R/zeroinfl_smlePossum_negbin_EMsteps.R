#' @importFrom bizicount zic.reg
E_step_zeroinfl_nb = function(prev_beta_fit = NULL, beta_cols, eta_cols, Y, ## parameters / variables for the outcome model Y|X,Z
                              beta = NULL, eta = NULL, theta = NULL, ## (optional) fixed parameters for outcome model Y|X,Z
                              prev_p, Bspline, comp_dat_unval, ## parameters / variables for the exposure model X|X*(,Z)
                              m, N, n, use_predict_pYgivX = TRUE) {
  ## sample sizes (for indexing)
  ## Update the phi_xi = P(X=x|Yi,Xi*,Z) for unvalidated subjects --------------
  ### Analysis model: P(Y|X,Z) -------------------------------------------------
  #### Calculate P(Y|X,Z) from negative binomial distribution ------------------
  #### Because of the latent zero-inflation component, use predicted probs -----
  if (use_predict_pYgivX) {
    pYgivX = predict(object = prev_beta_fit,
                     newdata = data.frame(comp_dat_unval),
                     type = "prob",
                     y.new = data.frame(comp_dat_unval)[, Y])
  } else {
    pYgivX = get_pYgivX_from_fit(fit = prev_beta_fit,
                                 beta = beta,
                                 eta = eta,
                                 theta = theta,
                                 data = comp_dat_unval,
                                 beta_cols = beta_cols,
                                 eta_cols = eta_cols,
                                 Y = Y)
  }
  ### Error mechanism: P(X|X*,Z) ---------------------------------------------
  pX = prev_p[rep(seq(1, m), each = (N - n)), ] *
    comp_dat_unval[, Bspline]
  ##############################################################################
  ## Estimate conditional expectations -----------------------------------------
  psi_num = c(pYgivX) * pX
  ### Update denominator -------------------------------------------------------
  #### Sum up all rows per id (e.g. sum over xk) -------------------------------
  psi_denom = rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
  #### Then sum over the sn splines --------------------------------------------
  psi_denom = rowSums(psi_denom)
  #### Avoid NaN resulting from dividing by 0 ----------------------------------
  psi_denom[psi_denom == 0] = 1
  ### And divide them! ---------------------------------------------------------
  psi_t = psi_num / psi_denom
  ### Update the w_kyi for unvalidated subjects --------------------------------
  ### by summing across the splines/ columns of psi_t --------------------------
  w_t = rowSums(psi_t)
  #### Add indicators for validated rows ---------------------------------------
  phi_aug = c(rep(x = 1, times = n), w_t)
  ##############################################################################
  ## Return vector of weights and matrix of psi_t ------------------------------
  return(list(phi_aug = phi_aug,
              psi_t = psi_t))
}

M_step_zeroinfl_nb = function(phi_aug, psi_t, ## weights and quantities from the E-step
                              re_analysis_formula, comp_dat_all, prev_beta, prev_eta, prev_theta, ## to update parameters for the outcome model Y|X,Z
                              prev_p, p_val_num, ## to update parameters for the exposure model X|X*(,Z)
                              m, N, n, ## sample sizes (for indexing)
                              tol, iterlim = 100) { ## criterion for convergence
  ## Update beta using weighted zero-inflated negative binomial regression -----
  new_fit = suppressWarnings(
    zic.reg(
      fmla = as.formula(re_analysis_formula),
      data = data.frame(cbind(comp_dat_all, phi_aug)),
      weights = phi_aug,
      dist = "nbinom",
      optimizer = "nlm",
      iterlim = iterlim,
      starts = c(prev_beta, prev_eta, 1 / prev_theta),
    )
  )
  new_beta = matrix(data = new_fit$coef[grepl(pattern = "ct_", x = names(new_fit$coef))],
                    ncol = 1)
  new_eta = matrix(data = new_fit$coef[grepl(pattern = "zi_", x = names(new_fit$coef))],
                   ncol = 1)
  new_theta = new_fit$coef["Theta"]
  ## Check for beta convergence ------------------------------------------------
  beta_conv = abs(new_beta - prev_beta) < tol
  eta_conv = abs(new_eta - prev_eta) < tol
  theta_conv = abs(new_theta - prev_theta) < tol
  ##############################################################################
  ## Update {p_kj} -------------------------------------------------------------
  ### Update numerators by summing u_t over i = 1, ..., N ----------------------
  new_p_num = p_val_num +
    rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
  new_p = t(t(new_p_num) / colSums(new_p_num))
  ### Check for convergence ----------------------------------------------------
  p_conv = abs(new_p - prev_p) < tol
  ##############################################################################
  # Check for global convergence -----------------------------------------------
  all_conv = c(beta_conv, eta_conv, theta_conv, p_conv)
  prop_conv = mean(all_conv)
  ##############################################################################
  ## Return updated parameters and proportion converged ------------------------
  return(list(new_beta_fit = new_fit,
              new_beta = new_beta,
              new_eta = new_eta,
              new_theta = new_theta,
              new_p = new_p,
              prop_conv = prop_conv))
}

M_step_zeroinfl_nb_ponly = function(psi_t, ## weights and quantities from the E-step
                                    prev_p, p_val_num, ## to update parameters for the exposure model X|X*(,Z)
                                    m, N, n, ## sample sizes (for indexing)
                                    tol) { ## criterion for convergence
  ##############################################################################
  ## Update {p_kj} -------------------------------------------------------------
  ### Update numerators by summing u_t over i = 1, ..., N ----------------------
  new_p_num = p_val_num +
    rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
  new_p = t(t(new_p_num) / colSums(new_p_num))
  ### Check for convergence ----------------------------------------------------
  p_conv = abs(new_p - prev_p) < tol
  prop_conv = mean(p_conv)
  ##############################################################################
  ## Return updated parameters and proportion converged ------------------------
  return(list(new_p = new_p,
              prop_conv = prop_conv))
}
