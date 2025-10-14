#' Sieve maximum likelihood estimation for negative binomial regression problems with covariate measurement error
#' This function returns the sieve maximum likelihood estimates (SMLEs) for the negative binomial regression model with covariate measurement error from Lotspeich et al. (2025+)
#'
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be the negative binomial model outcome.
#' @param error_formula formula, covariate error model formula (or coercible to formula), a formula expression as for other regression models. The response should be the error-free version of the error-prone of the covariate, and the covariate should be the names of the B-spline columns.
#' @param data dataset containing at least the variables included in \code{error_formula} and \code{analysis_formula}.
#' @param no_se Indicator for whether standard errors are desired. Defaults to \code{no_se = FALSE}.
#' @param pert_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is \code{pert_scale = 1}.
#' @param tol Tolerance between iterations in the EM algorithm used to define convergence.
#' @param max_iter Maximum number of iterations allowed in the EM algorithm.
#' @param output character, level of fitted model output to be returned. Defaults to \code{output = "coeff"}, but \code{output = "all"} is also possible.
#' @return
#' \item{coefficients}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{vcov}{variance-covariance matrix for \code{coefficients} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' @export
#' @importFrom MASS glm.nb
smlePossum_nb = function(analysis_formula, error_formula, data, no_se = TRUE, pert_scale = 1,
                         tol = 1E-4, max_iter = 1000, output = "coeff") {
  ##############################################################################
  # Extract variable names from user-specified formulas ------------------------
  ## Transform to formulas -----------------------------------------------------
  analysis_formula = as.formula(analysis_formula)
  error_formula = as.formula(error_formula)
  ## Outcome model -------------------------------------------------------------
  Y = as.character(analysis_formula)[2] ## outcome
  X_val = as.character(error_formula)[2] ## error-free covariate
  C = setdiff(x = unlist(strsplit(x = gsub(pattern = " ",
                                           replacement = "",
                                           x = as.character(as.formula(analysis_formula))[3]),
                                  split = "+",
                                  fixed = TRUE)),
              y = X_val)
  ## Error model ---------------------------------------------------------------
  Bspline = unlist(strsplit(x = gsub(pattern = " ",
                                     replacement = "",
                                     x = as.character(as.formula(error_formula))[3]),
                            split = "+",
                            fixed = TRUE))
  ## Extract analysis model matrix (complete cases) for structure --------------
  analysis_mat = model.matrix(object = analysis_formula,
                              data = data)
  ### Rewrite model formulas using column names from the model matrices --------
  re_analysis_formula = paste0(Y, "~",
                               paste(colnames(analysis_mat)[-1],
                                     collapse = "+"))
  ##############################################################################
  # Prepare for algorithm ------------------------------------------------------
  ## Define validation indicator and sample sizes ------------------------------
  data[, "Validated"] = as.numeric(!is.na(data[, X_val])) ## validation indicator
  N = nrow(data) ## total sample size (Phase I)
  n = sum(data[, "Validated"]) ## validation study sample size (Phase II)
  ## Create row numbers (for user-supplied data) -------------------------------
  data[, "row_num"] = 1:N
  ## Reorder so that the n validated subjects are first ------------------------
  data = data[order(as.numeric(data[, "Validated"]), decreasing = TRUE), ]
  ## Check for B-spline error --------------------------------------------------
  sn = ncol(data[, Bspline])
  if(0 %in% colSums(data[c(1:n), Bspline])) {
    warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)
    if(output == "coeff") {
      return(list(model_coeff = data.frame(coeff = NA, se = NA),
                  vcov = NA,
                  converged = NA,
                  se_converged = NA,
                  converged_msg = "B-spline error"))
    } else {
      return(list(model_coeff = data.frame(coeff = NA, se = NA),
                  vcov = NA,
                  predicted = NA,
                  converged = NA,
                  se_converged = NA,
                  converged_msg = "B-spline error"))
    }
  }
  ##############################################################################
  # Create complete data -------------------------------------------------------
  ## Save distinct X from validation study -------------------------------------
  x_obs = data.frame(unique(data[1:n, c(X_val)]))
  x_obs = data.frame(x_obs[order(x_obs[, 1]), ])
  m = nrow(x_obs)
  x_obs_stacked = do.call(what = rbind,
                          args = replicate(n = (N - n),
                                           expr = x_obs,
                                           simplify = FALSE))
  x_obs_stacked = data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
  colnames(x_obs) = colnames(x_obs_stacked) = c(X_val)
  ## Save static (X*,X,Y,C) since they don't change ----------------------------
  comp_dat_val = data[c(1:n), c(Y, X_val, C, Bspline)]
  comp_dat_val = merge(x = comp_dat_val,
                       y = data.frame(x_obs, k = 1:m),
                       all.x = TRUE)
  comp_dat_val = cbind(comp_dat_val, int = 1) ### add intercept column
  comp_dat_val = comp_dat_val[, c("int", Y, X_val, C, Bspline, "k")]
  comp_dat_val = data.matrix(comp_dat_val)

  # (m x n)xd vectors of each (one column per person, one row per x) -----------
  comp_dat_unval = suppressWarnings(
    data.matrix(
      cbind(data[-c(1:n), c(Y, C, Bspline)],
            x_obs_stacked,
            k = rep(seq(1, m), each = (N - n)))
    )
  )
  comp_dat_unval = cbind(comp_dat_unval, int = 1) ### add intercept column
  comp_dat_unval = comp_dat_unval[, c("int", Y, X_val, C, Bspline, "k")]
  comp_dat_all = rbind(comp_dat_val, comp_dat_unval)
  ##############################################################################
  # Initialize analysis model parameters (beta) --------------------------------
  ## Set initial values for beta and theta -------------------------------------
  beta_cols = c("int", X_val, C)
  cc_fit = glm.nb(formula = re_analysis_formula,
                  data = data.frame(comp_dat_val))
  prev_beta = beta0 = matrix(data = cc_fit$coefficients,
                             ncol = 1)
  prev_theta = theta0 = cc_fit$theta
  # prev_beta = beta0 = matrix(data = 0,
  #                            nrow = length(beta_cols),
  #                            ncol = 1)
  # prev_theta = 1E-3
  ### Set initial values for B-spline coefficients {p_kj} ----------------------
  p_val_num = rowsum(x = comp_dat_val[, Bspline],
                     group = comp_dat_val[, "k"],
                     reorder = TRUE)
  prev_p = p0 =  t(t(p_val_num) / colSums(p_val_num))
  ##############################################################################
  # Estimate beta, k, and p_kj and eta using EM algorithm ----------------------
  ## Set parameters for algorithm convergence ----------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  ## Otherwise, begin EM algorithm ---------------------------------------------
  while(it <= max_iter & !CONVERGED) {
    ############################################################################
    # E Step -------------------------------------------------------------------
    E_step_res = E_step_nb(prev_beta = prev_beta,
                           prev_theta = prev_theta,
                           Y = Y,
                           beta_cols = beta_cols,
                           prev_p = prev_p,
                           Bspline = Bspline,
                           comp_dat_unval = comp_dat_unval,
                           m = m,
                           N = N,
                           n = n)
    ############################################################################
    # M Step -------------------------------------------------------------------
    M_step_res = M_step_nb(phi_aug = E_step_res$phi_aug,
                           psi_t = E_step_res$psi_t,
                           re_analysis_formula = re_analysis_formula,
                           comp_dat_all = comp_dat_all,
                           prev_beta = prev_beta,
                           prev_theta = prev_theta,
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
    prev_beta = M_step_res$new_beta
    prev_theta = M_step_res$new_theta
    prev_p = M_step_res$new_p
  }
  ### Name rows of coefficients before preparing to return (below)
  new_coeff = c(M_step_res$new_beta, M_step_res$new_theta)
  names(new_coeff) = c(beta_cols, "Dispersion")

  # ----------------------------------- Estimate beta and eta using EM algorithm
  ##############################################################################
  # Check convergence statuses -------------------------------------------------
  if(!CONVERGED) {
    ## Return final estimates and convergence information ----------------------
    coeff_df = data.frame(coeff = NA,
                          se = NA,
                          z = NA,
                          p = NA)
    colnames(coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    return(list(coefficients = coeff_df,
                vcov = NA,
                converged = FALSE,
                se_converged = NA,
                converged_msg = "max_iter reached"))

  } else {
    CONVERGED_MSG = "Converged"
  }
  # ------------------------------------------------- Check convergence statuses
  ##############################################################################
  # Create list and return results ---------------------------------------------
  ## Predict X | X*, C (if requested) ------------------------------------------
  if (output == "all") {
    ## Create matrix with columns: (x_j) x (p_kj)
    xj_wide = matrix(data = unlist(x_obs),
                     nrow = nrow(M_step_res$new_p),
                     ncol = ncol(M_step_res$new_p),
                     byrow = FALSE)
    xj_phat = xj_wide * M_step_res$new_p

    ## Calculate predicted X given error-prone X* and Z
    xhat = data[, X_val] ### initialize with validated X (when non-missing)
    for (i in which(is.na(xhat))) {
      xhat[i] = smle_predict_x(row_data = data[i, ],
                               bspline_coeff = xj_phat)
    }
  }
  ## Calculate SE (if requested) and return ------------------------------------
  if(no_se){
    ## Build table of outcome model coefficients -------------------------------
    coeff_df = data.frame(coeff = new_coeff,
                          se = NA,
                          z = NA,
                          p = NA)
    colnames(coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    ## Return, depending on output preference ----------------------------------
    if(output == "coeff") {
      ## Return coefficients ---------------------------------------------------
      return(list(coefficients = coeff_df,
                  vcov = matrix(data = NA,
                                nrow = nrow(coeff_df),
                                ncol = nrow(coeff_df)),
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG))
    } else {
      ## Calculate l(beta, theta, p) -------------------------------------------
      od_loglik_conv =  smlePossum_negbin_od_ll(N = N,
                                                n = n,
                                                Y = Y,
                                                beta_cols = beta_cols,
                                                Bspline = Bspline,
                                                comp_dat_all = comp_dat_all,
                                                beta = M_step_res$new_beta,
                                                theta = M_step_res$new_theta,
                                                p = M_step_res$new_p)
      ## Return coefficients, B-spline coefficients, predictions ---------------
      return(list(coefficients = coeff_df,
                  bspline_coefficients = cbind(x_obs, M_step_res$new_p),
                  vcov = matrix(data = NA,
                                nrow = length(new_coeff),
                                ncol = length(new_coeff)),
                  predicted = xhat[order(data$orig_row)],
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG,
                  iterations = it,
                  od_loglik_at_conv = od_loglik_conv))
    }
  } else {
    # Estimate Cov(theta) using profile likelihood -----------------------------
    h_N = pert_scale * N ^ ( - 1 / 2) # perturbation ---------------------------
    ## Calculate l(beta, theta, p) ---------------------------------------------
    od_loglik_conv = smlePossum_negbin_od_ll(N = N,
                                             n = n,
                                             Y = Y,
                                             beta_cols = beta_cols,
                                             Bspline = Bspline,
                                             comp_dat_all = comp_dat_all,
                                             beta = M_step_res$new_beta,
                                             theta = M_step_res$new_theta,
                                             p = M_step_res$new_p)
    ## Setup information matrix with l(beta, theta) ----------------------------
    I_coeff = matrix(data = od_loglik_conv,
                     nrow = length(new_coeff),
                     ncol = length(new_coeff))
    ## Calculate single perturbation profile likelihoods pl(beta, theta) -------
    single_pert_theta = sapply(X = seq(1, ncol(I_coeff)),
                               FUN = smlePossum_negbin_pl,
                               beta_theta = new_coeff,
                               h_N = h_N,
                               N = N,
                               n = n,
                               Y = Y,
                               beta_cols = beta_cols,
                               Bspline = Bspline,
                               comp_dat_all = comp_dat_all,
                               p0 = M_step_res$new_p,
                               p_val_num = p_val_num,
                               tol = tol,
                               max_iter = max_iter)

    if (any(is.na(single_pert_theta))) {
      I_coeff = matrix(data = NA,
                       nrow = nrow(new_coeff),
                       ncol = nrow(new_coeff))
      SE_CONVERGED = FALSE
    } else {
      spt_wide = matrix(data = rep(c(single_pert_theta),
                                   times = ncol(I_coeff)),
                        ncol = ncol(I_coeff),
                        byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_coeff
      I_coeff = I_coeff - spt_wide - t(spt_wide)
      SE_CONVERGED = TRUE
    }

    for (c in 1:ncol(I_coeff)) {
      pert_coeff = new_coeff
      pert_coeff[c] = pert_coeff[c] + h_N
      double_pert_coeff = sapply(X = seq(c, ncol(I_coeff)),
                                 FUN = smlePossum_negbin_pl,
                                 beta_theta = pert_coeff,
                                 h_N = h_N,
                                 N = N,
                                 n = n,
                                 Y = Y,
                                 beta_cols = beta_cols,
                                 Bspline = Bspline,
                                 comp_dat_all = comp_dat_all,
                                 p0 = M_step_res$new_p,
                                 p_val_num = p_val_num,
                                 max_iter = max_iter,
                                 tol = tol)
      dpt = matrix(data = 0,
                   nrow = nrow(I_coeff),
                   ncol = ncol(I_coeff))
      dpt[c,c] = double_pert_coeff[1] # Put double on the diagonal
      if(c < ncol(I_coeff)) {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] = dpt[-(1:c), c] = double_pert_coeff[-1]
      }
      I_coeff = I_coeff + dpt
    }

    I_coeff = h_N ^ (- 2) * I_coeff

    cov_coeff = tryCatch(expr = - solve(I_coeff),
                         error = function(err) {
                           matrix(data = NA,
                                  nrow = nrow(I_coeff),
                                  ncol = ncol(I_coeff))
                         }
    )
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # if(any(diag(cov_coeff) < 0)) {
    #   warning("Negative variance estimate. Increase the pert_scale parameter and repeat variance estimation.")
    #   SE_CONVERGED = FALSE
    # }
    se_coeff = tryCatch(expr = sqrt(diag(cov_coeff)),
                        warning = function(w) {
                          matrix(NA, nrow = nrow(cov_coeff))
                        }
    )
    if (any(is.na(se_coeff))) { SE_CONVERGED = FALSE} else { TRUE }

    ## Return final estimates and convergence information ----------------------
    coeff_df = data.frame(coeff = new_coeff,
                          se = se_coeff,
                          z = new_coeff / se_coeff,
                          p = 2 * pnorm(q = abs(new_coeff / se_coeff), lower.tail = FALSE))
    rownames(coeff_df) = c(beta_cols, "Dispersion")
    colnames(coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(cov_coeff) = colnames(cov_coeff) = c(beta_cols, "Dispersion")

    if(output == "coeff") {
      return(list(coefficients = coeff_df,
                  vcov = cov_coeff,
                  converged = CONVERGED,
                  se_converged = SE_CONVERGED,
                  converged_msg = CONVERGED_MSG))
    } else {
      return(list(coefficients = coeff_df,
                  bspline_coefficients = cbind(x_obs, M_step_res$new_p),
                  vcov = cov_coeff,
                  predicted = xhat[order(data$orig_row)],
                  converged = CONVERGED,
                  se_converged = SE_CONVERGED,
                  converged_msg = CONVERGED_MSG,
                  iterations = it,
                  od_loglik_at_conv = od_loglik_conv))
    }
  }
  # --------------------------------------------- Create list and return results
}
