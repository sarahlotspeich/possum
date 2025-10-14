smlePossum_negbin_od_ll = function(N, n, Y, beta_cols, Bspline,
                                   comp_dat_all, beta, theta, p) {
  m = nrow(p)
  
  ##############################################################################
  # For validated subjects (first n rows) --------------------------------------
  ## Sum over log[P_theta(Yi|Xi)] ----------------------------------------------
  #### mu = exp(beta0 + beta1X + beta2Z + ) ...
  mu_beta = exp(as.numeric(comp_dat_all[c(1:n), beta_cols] %*% beta))
  #### Calculate P(Y|X,Z) from negative binomial distribution ------------------
  pYgivX = dnbinom(x = comp_dat_all[c(1:n), Y], 
                   size = theta, 
                   prob = (theta / (mu_beta + theta)))
  return_loglik = sum(log(pYgivX))
  ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ------------------------------------------
  pX = p[comp_dat_all[c(1:n), "k"], ]
  log_pX = log(pX)
  log_pX[log_pX == -Inf] = 0
  return_loglik = return_loglik + 
    sum(comp_dat_all[c(1:n), Bspline] * log_pX)
  
  ##############################################################################
  # For unvalidated subjects (last N - n rows) ---------------------------------
  ## Calculate P_beta,theta(y|x) for all (y,xk) --------------------------------
  #### mu = exp(beta0 + beta1X + beta2Z + ) ...
  mu_beta = exp(as.numeric(comp_dat_all[-c(1:n), beta_cols] %*% beta))
  #### Calculate P(Y|X,Z) from negative binomial distribution ------------------
  pYgivX = dnbinom(x = comp_dat_all[-c(1:n), Y], 
                   size = theta, 
                   prob = (theta / (mu_beta + theta)))
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
