#' Sieve maximum likelihood estimation for Poisson regression problems with covariate misclassification
#' This function returns the maximum likelihood estimates (MLEs) for the Poisson regression model with covariate misclassification from Mullan et al. (2024+)
#'
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be the Poisson model outcome, and, if needed, the offset can be provided as an \code{offset()} term.
#' @param family analysis model family, to be passed through to \code{glm}. See \code{?glm} for options.
#' @param error_formula formula, covariate error model formula (or coercible to formula), a formula expression as for other regression models. The response should be the error-free version of the error-prone of the covariate, and the covariate should be the names of the B-spline columns.
#' @param data dataset containing at least the variables included in \code{error_formula} and \code{analysis_formula}.
#' @param beta_init Initial values used to fit \code{analysis_formula}. Choices include (1) \code{"Zero"} (non-informative starting values, the default) or (2) \code{"Complete-data"} (estimated based on validated data only).
#' @param noSE Indicator for whether standard errors are desired. Defaults to \code{noSE = FALSE}.
#' @param hN_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is \code{hN_scale = 1}.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return
#' \item{coefficients}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{vcov}{variance-covariance matrix for \code{coefficients} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' @export
#' @importFrom numDeriv hessian
#' @importFrom MASS glm.nb
smlePossum_nb = function(analysis_formula, error_formula, data,
                         beta_init = "Zero", noSE = TRUE, hN_scale = 1, TOL = 1E-4, MAX_ITER = 1000) {
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
    if(output == "logORs") {
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
  ## If invalid beta_init specified, default to beta0 = "Zero" -----------------
  if(!(beta_init %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values for analysis model provided. Non-informative zeros assumed.")
    beta_init = "Zero"
  }
  ## Set initial values for beta -----------------------------------------------
  # #### Take some information from the complete-case fit (Poisson)
  # cc_fit = glm(formula = as.formula(re_analysis_formula),
  #              data = data.frame(comp_dat_val))
  # cc_fit = glm.nb(formula = as.formula(re_analysis_formula),
  #                 data = data.frame(comp_dat_val))
  # beta_cols = names(cc_fit$coefficients) ## column names
  # if(beta_init == "Complete-data") {
  #   prev_beta = beta0 = matrix(data = cc_fit$coefficients,
  #                              ncol = 1)
  # }
  beta_cols = c("int", X_val, C)
  if(beta_init == "Zero") {
    prev_beta = beta0 = matrix(data = 0,
                               nrow = length(beta_cols),
                               ncol = 1)
    prev_theta = 1E-3
  }
  ### Set initial values for B-spline coefficients {p_kj} ----------------------
  p_val_num = rowsum(x = comp_dat_val[, Bspline], 
                     group = comp_dat_val[, "k"], 
                     reorder = TRUE)
  prev_p = p0 =  t(t(p_val_num) / colSums(p_val_num))
  
  theta_design_mat = cbind(int = 1, 
                           comp_dat_all[, c(beta_cols)])
  ##############################################################################
  # Estimate beta, k, and p_kj and eta using EM algorithm ----------------------
  ## Set parameters for algorithm convergence ----------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  ## Otherwise, begin EM algorithm ---------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
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
    beta_conv = abs(new_beta - prev_beta) < TOL
    theta_conv = abs(new_theta - prev_theta) < TOL
    ############################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N -----------
    new_p_num = p_val_num +
      rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    new_p = t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv = abs(new_p - prev_p) < TOL
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
                converged_msg = "MAX_ITER reached"))

  } else {
    CONVERGED_MSG = "Converged"
  }
  # ------------------------------------------------- Check convergence statuses
  ##############################################################################
  # Create list and return results ---------------------------------------------
  if(noSE){
    ## Return final estimates and convergence information ----------------------
    coeff_df = data.frame(coeff = new_coeff,
                          se = NA,
                          z = NA,
                          p = NA)
    colnames(coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    return(list(coefficients = coeff_df,
                vcov = NA,
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG))
  } else {

    ## If alternative SE is preferred
    if(alternative_SE) {

      ### compute the Hessian
      hessian <- numDeriv::hessian(func = mle_loglik_nd,
                                   x = c(as.vector(new_beta), as.vector(new_eta)),
                                   method = "Richardson",
                                   beta_cols = beta_cols,
                                   eta_cols = eta_cols,
                                   Y = Y,
                                   offset = offset,
                                   X_unval = X_unval,
                                   X = X,
                                   comp_dat_val = comp_dat_val,
                                   comp_dat_unval = comp_dat_unval,
                                   noFN = noFN)

      ### use the Hessian to compute the standard error
      cov <- solve(hessian * - 1) #negate and invert Hessian for vcov @ MLE
      se <- sqrt(diag(cov)) #extract the standard errors

      SE_CONVERGED = !any(is.na(se))
      ### Split standard into the analysis and error model parameters ---------
      se_beta = se[c(1:nrow(prev_beta))]
      se_eta = se[-c(1:nrow(prev_beta))]
    } else {
    ## Calculate Cov(beta, eta) using numerical differentiation ----------------
    hN = hN_scale * N ^ ( - 1 / 2) # perturbation ------------------------------

    ## Calculate l(beta, eta) --------------------------------------------------
    od_loglik = mle_loglik(Y = Y,
                           offset = offset,
                           X_unval = X_unval,
                           X = X,
                           comp_dat_val = comp_dat_val,
                           comp_dat_unval = comp_dat_unval,
                           beta = new_beta,
                           eta = new_eta,
                           beta_cols = beta_cols,
                           eta_cols = eta_cols,
                           noFN = noFN)

    ## Calculate I(beta, eta) using numerical differentiation ------------------
    theta = rbind(new_beta, new_eta) ### create combined theta = (beta, eta)
    I = matrix(data = od_loglik,
               nrow = nrow(theta),
               ncol = nrow(theta))

    ### Compute log-likelihoods after single perturbations of beta and eta -----
    single_pert = vector(length = ncol(I))
    for (i in 1:length(single_pert)) {
      #### Define perturbed theta vector
      theta_pert = theta
      theta_pert[i] = theta_pert[i] + hN
      #### Calculate log-likelihood with perturbed theta vector
      od_loglik_pert = mle_loglik(Y = Y,
                                  offset = offset,
                                  X_unval = X_unval,
                                  X = X,
                                  comp_dat_val = comp_dat_val,
                                  comp_dat_unval = comp_dat_unval,
                                  beta = as.matrix(theta_pert[c(1:nrow(new_beta))]),
                                  eta = as.matrix(theta_pert[-c(1:nrow(new_beta))]),
                                  beta_cols = beta_cols,
                                  eta_cols = eta_cols,
                                  noFN = noFN)
      single_pert[i] = od_loglik_pert
    }

    ### Check for any elements that didn't converge ----------------------------
    if (any(is.na(single_pert))) {
      I = matrix(data = NA,
                 nrow = nrow(I),
                 ncol = nrow(I))
    } else {
      ### Create wide version of single perturbations --------------------------
      single_pert_wide = matrix(data = rep(x = single_pert,
                                           times = ncol(I)),
                                ncol = ncol(I),
                                byrow = FALSE)

      ### Using single_pert_wide, subtract the kth single perturbation from ----
      #### the kth row and kth column of the information matrix ----------------
      I = I - single_pert_wide - t(single_pert_wide)

      ### Compute log-likelihoods after double perturbations of beta and eta ---
      #### Create matrix of zeros for them (to be added to I_beta) -------------
      double_pert = matrix(data = 0,
                           nrow = nrow(I),
                           ncol = ncol(I))
      #### Loop over the rows and columns to fill double_pert ------------------
      for (r in 1:nrow(double_pert)) {
        #### First perturb the rth element in theta
        theta_pert = theta
        theta_pert[r] = theta_pert[r] + hN
        #### Then loop over further perturbing all elements in beta
        for (c in r:ncol(double_pert)) {
          #### Further perturb the cth element in theta
          theta_pert[c] = theta_pert[c] + hN
          #### Calculate log-likelihood with perturbed theta vector
          od_loglik_pert = mle_loglik(Y = Y,
                                      offset = offset,
                                      X_unval = X_unval,
                                      X = X,
                                      comp_dat_val = comp_dat_val,
                                      comp_dat_unval = comp_dat_unval,
                                      beta = as.matrix(theta_pert[c(1:nrow(new_beta))]),
                                      eta = as.matrix(theta_pert[-c(1:nrow(new_beta))]),
                                      beta_cols = beta_cols,
                                      eta_cols = eta_cols,
                                      noFN = noFN)
          double_pert[r, c] = double_pert[c, r] = od_loglik_pert
          #### Un-perturb the cth element in theta before incrementing c
          theta_pert[c] = theta_pert[c] - hN
        }
      }

      ### Add double perturbations to matrix of single perturbations -----------
      I = I + double_pert
    }

    ### Re-scale matrix of derivatives by the squared size of perturbation -----
    I = - hN ^ (- 2) * I

    cov = tryCatch(expr = solve(I),
                   error = function(err) {
                     matrix(data = NA,
                            nrow = nrow(I),
                            ncol = ncol(I))
                     }
                   )
    ## ---------------- Calculate Cov(beta, eta) using numerical differentiation
    ############################################################################
    ## Take square root of the diagonal elements to get standard errors --------
    se = tryCatch(expr = sqrt(diag(cov)),
                  warning = function(w) {
                    matrix(data = NA,
                           nrow = nrow(I))
                    }
                  )
    SE_CONVERGED = !any(is.na(se))
    ### Split them into the analysis and error model parameters ----------------
    se_beta = se[c(1:nrow(prev_beta))]
    se_eta = se[-c(1:nrow(prev_beta))]
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
