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
  prev_beta = beta0 = matrix(data = 0,
                             nrow = length(beta_cols),
                             ncol = 1)
  prev_theta = 1E-3
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
    # E Step -------------------------------------------------------------------
    ## Update the phi_xi = P(X=x|Yi,Xi*,Z) for unvalidated subjects ------------
    ### Analysis model: P(Y|X,Z) -----------------------------------------------
    #### mu = exp(beta0 + beta1X + beta2Z + ) ...
    mu_beta = exp(as.numeric(comp_dat_unval[, beta_cols] %*% prev_beta))
    #### Calculate P(Y|X,Z) from negative binomial distribution ----------------
    pYgivX = dnbinom(x = comp_dat_unval[, Y],
                     size = prev_theta,
                     prob = (prev_theta / (mu_beta + prev_theta)))
    ### Error mechanism: P(X|X*,Z) ---------------------------------------------
    pX = prev_p[rep(seq(1, m), each = (N - n)), ] *
      comp_dat_unval[, Bspline]
    ############################################################################
    ## Estimate conditional expectations ---------------------------------------
    psi_num = c(pYgivX) * pX
    ### Update denominator ------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk) ------------------
    psi_denom = rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
    #### Then sum over the sn splines -------------------------------
    psi_denom = rowSums(psi_denom)
    #### Avoid NaN resulting from dividing by 0 ---------------------
    psi_denom[psi_denom == 0] = 1
    ### And divide them! --------------------------------------------
    psi_t = psi_num / psi_denom
    ### Update the w_kyi for unvalidated subjects -------------------
    ### by summing across the splines/ columns of psi_t -------------
    w_t = rowSums(psi_t)
    #### Add indicators for validated rows -------------------------------------
    phi_aug = c(rep(x = 1, times = n), w_t)
    ############################################################################
    # M Step -------------------------------------------------------------------
    ## Update beta using weighted Poisson regression ---------------------------
    new_fit = glm.nb(formula = re_analysis_formula,
                     data = data.frame(cbind(comp_dat_all, phi_aug)),
                     weights = phi_aug)
    new_beta = matrix(data = new_fit$coefficients,
                      ncol = 1)
    new_theta = new_fit$theta
    ## Check for beta convergence ----------------------------------------------
    beta_conv = abs(new_beta - prev_beta) < tol
    theta_conv = abs(new_theta - prev_theta) < tol
    ############################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N -----------
    new_p_num = p_val_num +
      rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    new_p = t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv = abs(new_p - prev_p) < tol
    ############################################################################
    # Check for global convergence ---------------------------------------------
    all_conv = c(beta_conv, theta_conv, p_conv)
    CONVERGED = mean(all_conv) == 1
    # Update values for next iteration  ----------------------------------------
    it = it + 1
    prev_beta = new_beta
    prev_theta = new_theta
    prev_p = new_p
  }
  ### Name rows of coefficients before preparing to return (below)
  new_coeff = c(new_beta, new_theta)
  names(new_coeff) = c(beta_cols, "dispersion")

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
                     nrow = nrow(new_p),
                     ncol = ncol(new_p),
                     byrow = FALSE)
    xj_phat = xj_wide * new_p

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
                                ncol = nrow(coeff_df),
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG))
    } else {
      ## Calculate l(beta, theta, p) -------------------------------------------
      od_loglik_theta =  smlePossum_negbin_od_ll(N = N,
                                                 n = n,
                                                 Y = Y,
                                                 beta_cols = beta_cols,
                                                 Bspline = Bspline,
                                                 comp_dat_all = comp_dat_all,
                                                 beta = new_beta,
                                                 theta = new_theta,
                                                 p = new_p)
      ## Return coefficients, B-spline coefficients, predictions ---------------
      return(list(model_coeff = data.frame(coeff = new_theta,
                                           se = NA),
                  bspline_coeff = cbind(x_obs, new_p),
                  vcov = matrix(data = NA,
                                nrow = length(new_theta),
                                ncol = length(new_theta)),
                  predicted = xhat[order(data$orig_row)],
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG,
                  iterations = it,
                  od_loglik_at_conv = od_loglik_theta))
    }
  } else {
    # Estimate Cov(theta) using profile likelihood -----------------------------
    h_N = pert_scale * N ^ ( - 1 / 2) # perturbation ---------------------------
    ## Calculate l(beta, theta, p) ---------------------------------------------
    od_loglik_theta =  smlePossum_negbin_od_ll(N = N,
                                               n = n,
                                               Y = Y,
                                               beta_cols = beta_cols,
                                               Bspline = Bspline,
                                               comp_dat_all = comp_dat_all,
                                               beta = new_beta,
                                               theta = new_theta,
                                               p = new_p)
    ## Setup information matrix with l(beta, theta, p) -------------------------
    I_theta = matrix(data = od_loglik_theta,
                     nrow = (nrow(new_beta) + 1),
                     ncol = (nrow(new_beta) + 1))
    ## Calculate single perturbation profile likelihoods pl(beta, theta, p) ----
    single_pert_theta = sapply(X = seq(1, ncol(I_theta)),
                               FUN = pl_theta,
                               theta = new_theta,
                               h_N = h_N,
                               N = N,
                               n = n,
                               Y = Y,
                               X_val = X_val,
                               C = C,
                               Bspline = Bspline,
                               comp_dat_all = comp_dat_all,
                               p0 = new_p,
                               p_val_num = p_val_num,
                               tol = tol,
                               max_iter = max_iter)

    if (any(is.na(single_pert_theta))) {
      I_theta = matrix(data = NA,
                       nrow = nrow(new_theta),
                       ncol = nrow(new_theta))
      SE_CONVERGED = FALSE
    } else {
      spt_wide = matrix(data = rep(c(single_pert_theta),
                                   times = ncol(I_theta)),
                        ncol = ncol(I_theta),
                        byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
      I_theta = I_theta - spt_wide - t(spt_wide)
      SE_CONVERGED = TRUE
    }

    for (c in 1:ncol(I_theta)) {
      pert_theta = new_theta
      pert_theta[c] = pert_theta[c] + h_N
      double_pert_theta = sapply(X = seq(c, ncol(I_theta)),
                                 FUN = pl_theta,
                                 theta = pert_theta,
                                 h_N = h_N,
                                 N = N,
                                 n = n,
                                 Y = Y,
                                 X_val = X_val,
                                 C = C,
                                 Bspline = Bspline,
                                 comp_dat_all = comp_dat_all,
                                 p0 = new_p,
                                 p_val_num = p_val_num,
                                 max_iter = max_iter,
                                 tol = tol)
      dpt = matrix(data = 0,
                   nrow = nrow(I_theta),
                   ncol = ncol(I_theta))
      dpt[c,c] = double_pert_theta[1] #Put double on the diagonal
      if(c < ncol(I_theta)) {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] = dpt[-(1:c), c] = double_pert_theta[-1]
      }
      I_theta = I_theta + dpt
    }

    I_theta = h_N ^ (- 2) * I_theta

    cov_theta = tryCatch(expr = - solve(I_theta),
                         error = function(err) {
                           matrix(data = NA,
                                  nrow = nrow(I_theta),
                                  ncol = ncol(I_theta))
                         }
    )
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # if(any(diag(cov_theta) < 0)) {
    #   warning("Negative variance estimate. Increase the pert_scale parameter and repeat variance estimation.")
    #   SE_CONVERGED = FALSE
    # }
    se_theta = tryCatch(expr = sqrt(diag(cov_theta)),
                        warning = function(w) {
                          matrix(NA, nrow = nrow(prev_theta))
                        }
    )
    if (any(is.na(se_theta))) { SE_CONVERGED = FALSE} else { TRUE }
    if(output == "coeff") {
      return(list(model_coeff = data.frame(coeff = new_theta,
                                           se = se_theta),
                  vcov = cov_theta,
                  converged = CONVERGED,
                  se_converged = SE_CONVERGED,
                  converged_msg = CONVERGED_MSG))
    } else {
      return(list(model_coeff = data.frame(coeff = new_theta,
                                           se = se_theta),
                  bspline_coeff = cbind(x_obs, new_p),
                  vcov = cov_theta,
                  predicted = xhat[order(data$orig_row)],
                  converged = CONVERGED,
                  se_converged = SE_CONVERGED,
                  converged_msg = CONVERGED_MSG,
                  iterations = it,
                  od_loglik_at_conv = od_loglik_theta))
    }

    ## Return final estimates and convergence information ----------------------
    coeff_df = data.frame(coeff = new_beta,
                          se = se_beta,
                          z = new_beta / se_beta,
                          p = 2 * pnorm(q = abs(new_beta / se_beta), lower.tail = FALSE))
    rownames(coeff_df) = beta_cols
    colnames(coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    misclass_coeff_df = data.frame(coeff = new_eta,
                          se = se_eta,
                          z = new_eta / se_eta,
                          p = 2 * pnorm(q = abs(new_eta / se_eta), lower.tail = FALSE))
    rownames(misclass_coeff_df) = eta_cols
    colnames(misclass_coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    rownames(cov) = colnames(cov) = c(beta_cols, eta_cols)

    return(list(coefficients = coeff_df,
                misclass_coefficients = misclass_coeff_df,
                vcov = cov,
                converged = CONVERGED,
                se_converged = SE_CONVERGED,
                converged_msg = CONVERGED_MSG))
  }
  # --------------------------------------------- Create list and return results
}
