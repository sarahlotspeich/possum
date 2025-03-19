#' Maximum likelihood estimation for Poisson regression problems with covariate misclassification
#' This function returns the maximum likelihood estimates (MLEs) for the Poisson regression model with covariate misclassification from Mullan et al. (2024+)
#'
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be the Poisson model outcome, and, if needed, the offset can be provided as an \code{offset()} term.
#' @param family analysis model family, to be passed through to \code{glm}. See \code{?glm} for options.
#' @param error_formula misclassification model formula (or coercible to formula), a formula expression as for other regression models. The response should be the error-free version of the error-prone of the covariate, and the first predictor should be the error-prone version. 
#' @param data dataset containing at least the variables included in \code{error_formula} and \code{analysis_formula}.
#' @param beta_init Initial values used to fit \code{analysis_formula}. Choices include (1) \code{"Zero"} (non-informative starting values, the default) or (2) \code{"Complete-data"} (estimated based on validated data only).
#' @param eta_init Initial values used to fit \code{error_formula}. Choices include (1) \code{"Zero"} (non-informative starting values, the default) or (2) \code{"Complete-data"} (estimated based on validated data only).
#' @param noSE Indicator for whether standard errors are desired. Defaults to \code{noSE = FALSE}.
#' @param alternative_SE Indicator for whether alternative standard error should be calculated. Defaults to \code{alternative_SE = FALSE}
#' @param hN_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is \code{hN_scale = 1}.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return
#' \item{coefficients}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{misclass_coefficients}{dataframe with final coefficient and standard error estimates (where applicable) for the error model.}
#' \item{vcov}{variance-covariance matrix for \code{coefficients} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' @export
#' @importFrom numDeriv hessian
mlePossum = function(analysis_formula, family = poisson, error_formula, data,
                     beta_init = "Zero", eta_init = "Zero",
                     noSE = TRUE, alternative_SE = FALSE,
                     hN_scale = 1, TOL = 1E-4, MAX_ITER = 1000) {
  ## Extract variable names from user-specified formulas + model matrices ------
  ### Analysis model 
  Y = as.character(as.formula(analysis_formula))[2] ### Outcome 
  analysis_mat = model.matrix(object = analysis_formula, 
                              data = data) ### Design matrix for analysis formula
  ### Error model 
  X = as.character(as.formula(error_formula))[2] ### Validated covariate
  error_mat = model.matrix(object = error_formula, 
                           data = data)
  error_covar = colnames(error_mat)[-1] ### Exclude intercept 
  X_unval = error_covar[1] 
  
  ### Offset for the analysis model 
  offset = sub(pattern = "\\).*",
               replacement = "",
               x = sub(pattern = ".*offset\\(",
                       replacement = "",
                       x = grep(pattern = "offset", 
                                x = as.character(as.formula(analysis_formula)), 
                                value = TRUE)))
  offset = sub(pattern = "log\\(", replacement = "", x = offset) ### Check for log()
  if(length(offset) == 0) offset = NULL ## fixes data typing issue if no offset
  
  ### Rewrite model formulas using column names from the model matrices --------
  re_analysis_formula = paste0(Y, "~", 
                               paste(colnames(analysis_mat)[-1], 
                                     collapse = "+")) 
  if (!is.null(offset)) {
    re_analysis_formula = paste0(re_analysis_formula, "+", 
                                 sub(pattern = ".*offset\\(",
                                     replacement = "offset\\(",
                                     x = grep(pattern = "offset", 
                                              x = as.character(as.formula(analysis_formula)), 
                                              value = TRUE)))
  }
  re_error_formula = paste0(X, "~", 
                            paste(colnames(error_mat)[-1], 
                                  collapse = "+")) 
  
  # Prepare for algorithm ------------------------------------------------------
  data[, "Validated"] = as.numeric(!is.na(data[, X])) ## validation indicator
  N = nrow(data) ## total sample size (Phase I)
  n = sum(data[, "Validated"]) ## validation study sample size (Phase II)

  ## Reorder so that the n validated subjects are first ------------------------
  data = data[order(as.numeric(data[, "Validated"]), decreasing = TRUE), ]

  # Check for possibility of false negatives in validated data -----------------
  noFN = !any(data[1:n, X] == 1 & data[1:n, X_unval] == 0)

  ## Create row numbers --------------------------------------------------------
  data[, "row_num"] = 1:N
  ##############################################################################
  ## Save static (X*,X,Y,Z) for validated rows since they don't change ---------
  #comp_dat_val = data.matrix(data[c(1:n), c(Y, offset, X_unval, X, Z, "row_num")])
  comp_dat_val = make_complete_data(data = data, 
                                    analysis_formula = analysis_formula, 
                                    error_formula = error_formula, 
                                    rows = 1:n, 
                                    Y = Y, 
                                    X = X,
                                    offset = offset)

  ## Create augmented (X*,x,Y,Z) for unvalidated rows --------------------------
  ### First (N-n) rows assume X = 0 --------------------------------------------
  comp_dat0 = make_complete_data(data = data,
                                 analysis_formula = analysis_formula, 
                                 error_formula = error_formula, 
                                 rows = -c(1:n), 
                                 Y = Y, 
                                 X = X,
                                 offset = offset, 
                                 x = 0)
  # comp_dat0 = data.matrix(data[-c(1:n), c(Y, offset, X_unval, X, Z, "row_num")])
  # comp_dat0[, X] = 0
  ### Last (N-n) rows assume X = 1 ---------------------------------------------
  comp_dat1 = make_complete_data(data = data,
                                 analysis_formula = analysis_formula, 
                                 error_formula = error_formula, 
                                 rows = -c(1:n), 
                                 Y = Y, 
                                 X = X,
                                 offset = offset, 
                                 x = 1)
  # comp_dat1 = data.matrix(data[-c(1:n), c(Y, offset, X_unval, X, Z, "row_num")])
  # comp_dat1[, X] = 1
  ### Put them together --------------------------------------------------------
  comp_dat_unval = data.matrix(rbind(comp_dat0,
                                     comp_dat1))
  colnames(comp_dat_unval) = colnames(comp_dat_val) ## Coerce colnames to match
  
  ## Create augmented "complete" dataset of validated and unvalidated ----------
  comp_dat_all = data.matrix(rbind(comp_dat_val, comp_dat_unval))
  
  ##############################################################################
  ## Initialize analysis model parameters (beta) -------------------------------
  ### If invalid beta0 specified, default to beta0 = "Zero" --------------------
  if(!(beta_init %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values for analysis model provided. Non-informative zeros assumed.")
    beta_init = "Zero"
  }
  ### Set initial values for beta ----------------------------------------------
  if(beta_init == "Complete-data") {
    cc_fit = glm(formula = as.formula(re_analysis_formula),
                 family = family,
                 data = comp_dat_val)
    prev_beta = beta0 = matrix(data = cc_fit$coefficients,
                               ncol = 1)
    beta_cols = names(cc_fit$coefficients) ## column names
  } 
  if(beta_init == "Zero") {
    prev_beta = beta0 = matrix(data = 0,
                               nrow = nrow(prev_beta),
                               ncol = 1)
  }  
  ### Set initial values for eta -----------------------------------------------
  #### If one-sided (no false negatives), X_unval is not included in this model
  if (noFN) {
    subset_X_unval_one = comp_dat_val[comp_dat_val[, X_unval] == 1, ]
    message("Error model was modified to exclude unvalidated covariate, since errors are one-sided (false positives only).")
    if (sum(!grepl(pattern = X_unval, x = error_covar)) == 0) {
      re_error_formula = paste(X, "~ 1")
    } else {
      re_error_formula = paste(X, "~", 
                               paste(error_covar[!grepl(pattern = X_unval, x = error_covar)], 
                                     collapse = " + "))
    }
  } 
  if(eta_init == "Complete-data") {
    if (noFN) {
      cc_fit = glm(formula = as.formula(re_error_formula),
                   family = "binomial",
                   data = subset_X_unval_one)
    } else {
      cc_fit = glm(formula = as.formula(re_error_formula),
                   family = "binomial",
                   data = comp_dat_val)
    }
    prev_eta = eta0 = matrix(cc_fit$coefficients,
                             ncol = 1)
    eta_cols = names(cc_fit$coefficients) ## column names
  } 
  if(eta_init == "Zero") {
    prev_eta = eta0 = matrix(data = 0,
                             nrow = nrow(prev_eta),
                             ncol = 1)
  }
  # ------------------------------------------------------ Prepare for algorithm
  ##############################################################################
  # Estimate beta and eta using EM algorithm -----------------------------------
  ## Set parameters for algorithm convergence ----------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1

  # Check if we can even do algorithm ------------------------------------------
  queried_ppv = sum(data[,X] == 1 & data[,X_unval] == 1, na.rm = TRUE) /
    sum(data[,X_unval] == 1 & !is.na(data[,X]), na.rm = TRUE)

  ## If PPV among queried subset is almost perfect, just fit usual model -------
  if (round(queried_ppv, 3) == 1) {
    # vanilla_mod <- glm(formula = as.formula(analysis_formula),
    #                    family = family,
    #                    data = data)
    return(list(coefficients = data.frame(coeff = NA,
                                          se = NA),
                misclass_coefficients = data.frame(coeff = NA, se = NA),
                vcov = vcov(vanilla_mod),
                converged = NA,
                se_converged = NA,
                converged_msg = "Validated PPV = 1, use standard GLM "))
  }

  ## Otherwise, begin EM algorithm ---------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step -------------------------------------------------------------------
    ## Update the phi_xi = P(X=x|Yi,Xi*,Z) for unvalidated subjects ------------
    ### Analysis model: P(Y|X,Z) -----------------------------------------------
    #### mu = beta0 + beta1X + beta2Z + ...
    mu_beta = as.numeric(comp_dat_unval[, beta_cols] %*% prev_beta)
    #### lambda = exp(beta0 + beta1X + beta2Z + ... )
    lambda = exp(mu_beta)
    #### If offset specified, lambda = offset x exp(beta0 + beta1X + beta2Z + ... )
    if (!is.null(offset)) {
      lambda = comp_dat_unval[, offset] * lambda
    }
    #### Calculate P(Y|X) from Poisson distribution ----------------------------
    pYgivX = dpois(x = comp_dat_unval[, Y],
                   lambda = lambda)
    ############################################################################
    ### Misclassification mechanism: P(X|X*,Z) ---------------------------------
    #### mu = eta0 + eta1X* + eta2Z + ...
    mu_eta = as.numeric(comp_dat_unval[, eta_cols] %*% prev_eta)
    #### Calculate P(X|X*,Z) from Bernoulli distribution ---------------------
    pXgivXstar = dbinom(x = comp_dat_unval[, X],
                        size = 1,
                        prob = 1 / (1 + exp(- mu_eta)))
    #### Save min/max P(X|X,Z) to check for numerical 0/1 later --------------
    min_pXgivXstar = min(pXgivXstar)
    max_pXgivXstar = max(pXgivXstar)
    if (noFN) { #### If one-sided errors, logistic regression on just X*=1 -----
      #### Force P(X=0|X*=0,Z)=1 and P(X=1|X*=0,Z)=0 for all Z -----------------
      pXgivXstar[which(comp_dat_unval[, X_unval] == 0 & comp_dat_unval[, X] == 0)] = 1
      pXgivXstar[which(comp_dat_unval[, X_unval] == 0 & comp_dat_unval[, X] == 1)] = 0
    }
    ############################################################################
    ## Estimate conditional expectations ---------------------------------------
    ### Update numerator -------------------------------------------------------
    #### P(Y|X,Z)P(X|X*,Z) -----------------------------------------------------
    phi_num = pYgivX * pXgivXstar ##### dim: 2(N - n) x 1
    phi_num_wide = matrix(data = phi_num,
                          nrow = (N - n),
                          ncol = 2,
                          byrow = FALSE)
    ### Update denominator -----------------------------------------------------
    #### P(Y|X=0,Z)P(X=0|X*) + P(Y|X=1,Z)P(X=1|X*) -----------------------------
    phi_denom = rowSums(phi_num_wide) ##### dim: (N - n) x 1
    #### Avoid NaN resulting from dividing by 0 --------------------------------
    phi_denom[phi_denom == 0] = 1
    ### Divide them to get psi = E{I(X=x)|Y,X*} --------------------------------
    psi = phi_num / rep(x = phi_denom, times = 2)
    #### Add indicators for validated rows -------------------------------------
    phi_aug = c(rep(x = 1, times = n), psi)
    ############################################################################
    # M Step -------------------------------------------------------------------
    ## Update beta using weighted Poisson regression ---------------------------
    new_beta = suppressWarnings(
      matrix(data = glm(formula = re_analysis_formula,
                        family = family,
                        data = data.frame(cbind(comp_dat_all, phi_aug)),
                        weights = phi_aug)$coefficients,
             ncol = 1)
    )
    ## Check for beta convergence ----------------------------------------------
    beta_conv = abs(new_beta - prev_beta) < TOL
    ############################################################################
    ## Update eta using weighted logistic regression ---------------------------
    if (noFN) {
      which_unval_case = which(comp_dat_all[, X_unval] == 1)
      new_eta = suppressWarnings(
        matrix(data = glm(formula = re_error_formula,
                          family = binomial,
                          data = data.frame(cbind(comp_dat_all, phi_aug))[which_unval_case, ],
                          weights = phi_aug)$coefficients,
               ncol = 1)
      )
    } else {
      new_eta = suppressWarnings(
        matrix(data = glm(formula = re_error_formula,
                          family = binomial,
                          data = data.frame(comp_dat_all),
                          weights = phi_aug)$coefficients,
               ncol = 1)
      )
    }
    ## Check for beta convergence ----------------------------------------------
    eta_conv = abs(new_eta - prev_eta) < TOL
    ############################################################################
    # Check for global convergence ---------------------------------------------
    all_conv = c(beta_conv, eta_conv)
    CONVERGED = mean(all_conv) == 1
    # Update values for next iteration  ----------------------------------------
    it = it + 1
    prev_beta = new_beta
    prev_eta = new_eta
  }
  ### Name rows of coefficients before preparing to return (below) 
  beta_cols = beta_cols 
  rownames(new_eta) = eta_cols
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
    
    misclass_coeff_df = data.frame(coeff = NA,
                                   se = NA, 
                                   z = NA, 
                                   p = NA)
    colnames(misclass_coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    return(list(coefficients = coeff_df,
                misclass_coefficients = misclass_coeff_df,
                vcov = NA,
                converged = FALSE,
                se_converged = NA,
                converged_msg = "MAX_ITER reached"))
    
  } else {
    ## Even if algorithm converged, check for fitted probabilities close to ----
    ## Zero or one with the etas at convergence --------------------------------
    if (min_pXgivXstar < 1e-308 | max_pXgivXstar > (1-1e-16)) {
      CONVERGED_MSG = "Fitted probabilities numerically 0 or 1 at convergence"
    } else {
      CONVERGED_MSG = "Converged"
    }
  }
  # ------------------------------------------------- Check convergence statuses
  ##############################################################################
  # Create list and return results ---------------------------------------------
  if(noSE){
    ## Return final estimates and convergence information ----------------------
    coeff_df = data.frame(coeff = new_beta,
                          se = NA, 
                          z = NA, 
                          p = NA)
    colnames(coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    misclass_coeff_df = data.frame(coeff = new_eta,
                                   se = NA, 
                                   z = NA, 
                                   p = NA)
    colnames(misclass_coeff_df) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    return(list(coefficients = coeff_df,
                misclass_coefficients = misclass_coeff_df,
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
