#' @importFrom bizicount dzinb
zeroinfl_smlePossum_negbin_od_ll = function(N, n, prev_beta_fit = NULL, Y, Bspline, comp_dat_all,
                                            beta = NULL, eta = NULL, theta = NULL, beta_cols = NULL, eta_cols = NULL,
                                            p, use_predict_pYgivX = TRUE) {
  m = nrow(p)

  ##############################################################################
  # For validated subjects (first n rows) --------------------------------------
  ## Sum over log[P_theta(Yi|Xi)] ----------------------------------------------
  ### Calculate P(Y|X,Z) from zero-inflated negative binomial distribution -----
  ### Because of the latent zero-inflation component, use predicted probs ------
  if (use_predict_pYgivX) {
    pYgivX = predict(object = prev_beta_fit,
                     newdata = data.frame(comp_dat_all[c(1:n),]),
                     type = "prob",
                     y.new = data.frame(comp_dat_all[c(1:n), Y]))
  } else {
    pYgivX = get_pYgivX_from_fit(fit = prev_beta_fit,
                                 beta = beta,
                                 eta = eta,
                                 theta = theta,
                                 data = data.frame(comp_dat_all[c(1:n), ]),
                                 beta_cols = beta_cols,
                                 eta_cols = eta_cols,
                                 Y = Y)
  }
  return_loglik = sum(log(pYgivX))
  ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ------------------------------------------
  pX = p[comp_dat_all[c(1:n), "k"], ]
  log_pX = log(pX)
  log_pX[log_pX == -Inf] = 0
  return_loglik = return_loglik +
    sum(comp_dat_all[c(1:n), Bspline] * log_pX)
  ##############################################################################
  # For unvalidated subjects (last N - n rows) ---------------------------------
  ## Calculate P(Y|X,Z) from zero-inflated negative binomial distribution -----
  ### Because of the latent zero-inflation component, use predicted probs ------
  if (use_predict_pYgivX) {
    pYgivX = predict(object = prev_beta_fit,
                     newdata = data.frame(comp_dat_all[-c(1:n),]),
                     type = "prob",
                     y.new = data.frame(comp_dat_all[-c(1:n), Y]))
  } else {
    pYgivX = get_pYgivX_from_fit(fit = prev_beta_fit,
                                 beta = beta,
                                 eta = eta,
                                 theta = theta,
                                 data = data.frame(comp_dat_all[-c(1:n), ]),
                                 beta_cols = beta_cols,
                                 eta_cols = eta_cols,
                                 Y = Y)
  }
  ## Calculate Bj(Xi*) p_kj for all (k,j) --------------------------------------
  pX = rowSums(p[comp_dat_all[-c(1:n), "k"], ] * comp_dat_all[-c(1:n), Bspline])
  ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  ################################################################################
  ## Calculate sum of P(y|xk) x Bj(X*) x p_kj ------------------------------------
  person_sum = rowsum(pYgivX * pX,
                       group = rep(seq(1, (N - n)),
                                   times = m))
  log_person_sum = log(person_sum)
  log_person_sum[log_person_sum == -Inf] = 0
  ## And sum over them all -------------------------------------------------------
  return_loglik = return_loglik + sum(log_person_sum)
  return(return_loglik)
}
