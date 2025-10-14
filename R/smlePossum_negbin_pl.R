smlePossum_negbin_pl = function(k, beta_theta, h_N, N, n, Y, beta_cols, Bspline, comp_dat_all,
                                p0 = NULL, p_val_num = NULL, tol, max_iter) {
  # Perturb the kth entry in beta_theta (combined) by h_N
  pert = beta_theta
  pert[k] = pert[k] + h_N

  # Split perturbed vector into beta and theta
  pert_theta = pert[length(pert)]
  pert_beta = pert[-length(pert)]

  # Find the B-spline coefficients p based on pert
  pl_params = smlePossum_negbin_profile_out(beta = pert_beta,
                                            theta = pert_theta,
                                            N = N,
                                            n = n,
                                            Y = Y,
                                            beta_cols = beta_cols,
                                            Bspline = Bspline,
                                            comp_dat_all = comp_dat_all,
                                            p0 = p0,
                                            p_val_num = p_val_num,
                                            tol = tol,
                                            max_iter = max_iter)

  # If profile parameters converged, calculate observed-data log-likelihood
  if(pl_params$converged) {
    od_loglik_pert = smlePossum_negbin_od_ll(N = N,
                                             n = n,
                                             Y = Y,
                                             beta_cols = beta_cols,
                                             Bspline = Bspline,
                                             comp_dat_all = comp_dat_all,
                                             beta = pert_beta,
                                             theta = pert_theta,
                                             p = pl_params$p_at_conv)

  } else { od_loglik_pert = NA }
  return(od_loglik_pert)
}
