#' Maximum likelihood estimation for Poisson regression problems with covariate misclassification
#' This function returns the maximum likelihood estimates (MLEs) for the Poisson regression model with covariate misclassification from Mullan et al. (2024+)
#'
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be the Poisson model outcome, and, if needed, the offset can be provided as an \code{offset()} term.
#' @param error_formula misclassification model formula (or coercible to formula), a formula expression as for other regression models. The response should be the error-free version of the error-prone of the covariate.
#' @param data dataset containing at least the variables included in \code{error_formula} and \code{analysis_formula}.
#' @param beta_init Initial values used to fit \code{analysis_formula}. Choices include (1) \code{"Zero"} (non-informative starting values, the default) or (2) \code{"Complete-data"} (estimated based on validated data only).
#' @param eta_init Initial values used to fit \code{error_formula}. Choices include (1) \code{"Zero"} (non-informative starting values, the default) or (2) \code{"Complete-data"} (estimated based on validated data only).
#' @param noFN logical, if \code{noFN = FALSE} (the default), then it is assumed that there can be both false positives and false negatives in the error-prone exposure. If \code{noFN = TRUE}, the error mechanism is restricted to only false positives.
#' @param noSE Indicator for whether standard errors are desired. Defaults to \code{noSE = FALSE}.
#' @param hN_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is \code{hN_scale = 1}.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return 
#' \item{coefficients}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{misclass_coefficients}{dataframe with final coefficient and standard error estimates (where applicable) for the error model.}
#' \item{vcov}{variance-covariance matrix for \code{coeff} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' @export

mlePossum2 = function(analysis_formula, error_formula, data, beta_init = "Zero", eta_init = "Zero", noFN = FALSE, noSE = TRUE, hN_scale = 1, TOL = 1E-4, MAX_ITER = 1000) {
  ## Extract variable names from user-specified formulas
  Y = as.character(as.formula(analysis_formula))[2]
  X = as.character(as.formula(error_formula))[2]
  analysis_covar = unlist(strsplit(x = gsub(pattern = " ",
                                            replacement = "",
                                            x = as.character(as.formula(analysis_formula))[3]),
                                   split = "+",
                                   fixed = TRUE))
  error_covar = unlist(strsplit(x = gsub(pattern = " ",
                                         replacement = "",
                                         x = as.character(as.formula(error_formula))[3]),
                                split = "+",
                                fixed = TRUE))
  X_unval = setdiff(error_covar, analysis_covar)
  Z = intersect(error_covar, analysis_covar)
  offset = sub(pattern = "\\).*", 
               replacement = "", 
               x = sub(pattern = ".*\\(", 
                       replacement = "", 
                       x = setdiff(analysis_covar, c(X, Z))))
  
  # Prepare for algorithm ------------------------------------------------------
  data[, "Validated"] = as.numeric(!is.na(data[, X])) ## validation indicator
  N = nrow(data) ## total sample size (Phase I)
  n = sum(data[, "Validated"]) ## validation study sample size (Phase II)
  
  ## Reorder so that the n validated subjects are first ------------------------
  data = data[order(as.numeric(data[, "Validated"]), decreasing = TRUE), ]
  
  ## Create row numbers --------------------------------------------------------
  data[, "row_num"] = 1:N
  ##############################################################################
  ## Save static (X*,X,Y,Z) for validated rows since they don't change ---------
  comp_dat_val = data.matrix(data[c(1:n), c(Y, offset, X_unval, X, Z, "row_num")])
  
  ## Create augmented (X*,x,Y,Z) for unvalidated rows --------------------------
  ### First (N-n) rows assume X = 0 --------------------------------------------
  comp_dat0 = data.matrix(data[-c(1:n), c(Y, offset, X_unval, X, Z, "row_num")])
  comp_dat0[, X] = 0 
  ### Last (N-n) rows assume X = 1 ---------------------------------------------
  comp_dat1 = data.matrix(data[-c(1:n), c(Y, offset, X_unval, X, Z, "row_num")])
  comp_dat1[, X] = 1
  ### Put them together --------------------------------------------------------
  comp_dat_unval = rbind(comp_dat0, 
                         comp_dat1)
  
  ## Create augmented "complete" dataset of validated and unvalidated ----------
  comp_dat_all = rbind(comp_dat_val, comp_dat_unval)
  ##############################################################################
  ## Initialize analysis model parameters (beta) -------------------------------
  ### If invalid beta0 specified, default to beta0 = "Zero" --------------------
  if(!(beta_init %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values for analysis model provided. Non-informative zeros assumed.")
    beta_init = "Zero"
  }
  ### Set initial values for beta ----------------------------------------------
  if(beta_init == "Zero") {
    prev_beta = beta0 = matrix(data = 0, 
                               nrow = (length(c(X, Z)) + 1), 
                               ncol = 1)
  } else if(beta_init == "Complete-data") {
    prev_beta = beta0 = matrix(glm(formula = as.formula(analysis_formula),
                                   family = "poisson",
                                   data = data[c(1:n), ])$coefficients,
                               ncol = 1)
  }
  ### Set initial values for eta -----------------------------------------------
  #### If one-sided (no false negatives), X_unval is not included in this model
  if (noFN) {
    subset_X_unval_one = data[data[, X_unval] == 1, ] 
    message("Error model was modified to exclude unvalidated covariate, since errors are one-sided (false positives only).")
    error_model_covar = setdiff(x = error_covar, y = X_unval)
    if (length(error_model_covar) == 0) {
      error_formula = paste(X, "~ 1")
    } else {
      error_formula = paste(X, "~", paste(Z, collapse = " + "))
    }
  }
  if(eta_init == "Zero") {
    if (noFN) {
      #### If one-sided (no false negatives), X is not included in this model --
      prev_eta = eta0 = matrix(data = 0, 
                               nrow = (length(Z) + 1), 
                               ncol = 1)
    } else {
      prev_eta = eta0 = matrix(data = 0, 
                               nrow = (length(c(X, Z)) + 1), 
                               ncol = 1)
    }
  } else if(eta_init == "Complete-data") {
    if (noFN) {
      prev_eta = eta0 = matrix(glm(formula = as.formula(error_formula),
                                   family = "binomial",
                                   data = data[c(1:n), ], 
                                   subset = data[c(1:n), X_unval] == 1)$coefficients,
                               ncol = 1)
    } else {
      prev_eta = eta0 = matrix(glm(formula = as.formula(error_formula),
                                   family = "binomial",
                                   data = data[c(1:n), ])$coefficients,
                               ncol = 1)
    }
  }
  # ------------------------------------------------------ Prepare for algorithm
  ##############################################################################
  # Estimate beta and ppv using EM algorithm -----------------------------------
  ## Set parameters for algorithm convergence ----------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  ## Begin algorithm -----------------------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step -------------------------------------------------------------------
    ## Update the psi_xi = P(X=x|Yi,Xi*,Z) for unvalidated subjects ------------
    ### Analysis model: P(Y|X,Z) -----------------------------------------------
    #### mu = beta0 + beta1X + beta2Z + ... 
    mu_beta = as.numeric(cbind(int = 1, comp_dat_unval[, c(X, Z)]) %*% prev_beta)
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
    if (noFN) { #### If one-sided errors, logistic regression on just X*=1 -----
      #### mu = eta0 + eta1Z + ... 
      mu_eta = as.numeric(cbind(int = 1, comp_dat_unval[, Z]) %*% prev_eta)
      #### Calculate P(X|X*=1,Z) from Bernoulli distribution -------------------
      pXgivXstar = dbinom(x = comp_dat_unval[, X], 
                          size = 1, 
                          prob = 1 / (1 + exp(- mu_eta)))
      #### Force P(X=0|X*=0,Z)=1 and P(X=1|X*=0,Z)=0 for all Z -----------------
      pXgivXstar[which(comp_dat_unval[, X_unval] == 0 & comp_dat_unval[, X] == 0)] = 1
      pXgivXstar[which(comp_dat_unval[, X_unval] == 0 & comp_dat_unval[, X] == 0)] = 0
    } else { #### If two-sided errors, logistic regression on all rows ---------
      #### mu = eta0 + eta1X* + eta2Z + ... 
      mu_eta = as.numeric(cbind(int = 1, comp_dat_unval[, c(X_unval, Z)]) %*% prev_eta)
      #### Calculate P(X|X*,Z) from Bernoulli distribution ---------------------
      pXgivXstar = dbinom(x = comp_dat_unval[, X], 
                          size = 1, 
                          prob = 1 / (1 + exp(- mu_eta)))
    }
    ############################################################################
    ## Estimate conditional expectations ---------------------------------------
    ### Update numerator -------------------------------------------------------
    #### P(Y|X,Z)P(X|X*,Z) -----------------------------------------------------
    psi_num = pYgivX * pXgivXstar ##### dim: 2(N - n) x 1
    psi_num_wide = matrix(data = psi_num, 
                          nrow = (N - n), 
                          ncol = 2, 
                          byrow = FALSE)
    ### Update denominator -----------------------------------------------------
    #### P(Y|X=0,Z)P(X=0|X*) + P(Y|X=1,Z)P(X=1|X*) -----------------------------
    psi_denom = rowSums(psi_num_wide) ##### dim: (N - n) x 1
    #### Avoid NaN resulting from dividing by 0 --------------------------------
    psi_denom[psi_denom == 0] = 1
    ### Divide them to get psi = E{I(X=x)|Y,X*} --------------------------------
    psi = psi_num / rep(x = psi_denom, times = 2) 
    #### Add indicators for validated rows -------------------------------------
    psi_aug = c(rep(x = 1, times = n), psi)
    ############################################################################
    # M Step -------------------------------------------------------------------
    ## Update beta using weighted Poisson regression ---------------------------
    new_beta = suppressWarnings(
      matrix(data = glm(formula = analysis_formula,
                        family = poisson,
                        data = data.frame(cbind(comp_dat_all, psi_aug)),
                        weights = psi_aug)$coefficients,
             ncol = 1)
    )
    ## Check for beta convergence ----------------------------------------------
    beta_conv = abs(new_beta - prev_beta) < TOL
    ############################################################################
    ## Update eta using weighted logistic regression ---------------------------
    if (noFN) {
      which_unval_case = which(comp_dat_all[, X_unval] == 1)
      new_eta = suppressWarnings(
        matrix(data = glm(formula = error_formula,
                          family = binomial,
                          data = data.frame(cbind(comp_dat_all, psi_aug))[which_unval_case, ],
                          weights = psi_aug)$coefficients,
               ncol = 1)
      )
    } else {
      new_eta = suppressWarnings(
        matrix(data = glm(formula = error_formula,
                          family = binomial,
                          data = data.frame(cbind(comp_dat_all, psi_aug)),
                          weights = psi_aug)$coefficients,
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
  rownames(new_beta) = c("Intercept", X, Z)
  if(noFN) {
    rownames(new_eta) = c("Intercept", Z)
  } else {
    rownames(new_eta) = c("Intercept", X_unval, Z)
  }
  
  if(!CONVERGED) {
    if(it > MAX_ITER) {
      CONVERGED_MSG = "MAX_ITER reached"
    }
    
    return(list(coefficients = data.frame(coeff = NA, se = NA),
                misclass_coefficients = data.frame(coeff = NA, se = NA),
                vcov = NA,
                converged = FALSE,
                se_converged = NA,
                converged_msg = "MAX_ITER reached"))
  }
  
  if(CONVERGED) { 
    CONVERGED_MSG = "Converged" 
  }
  
  # ---------------------------------------------- Estimate beta using EM
  if(noSE){
    return(list(coefficients = data.frame(coeff = new_beta, se = NA),
                misclass_coefficients = data.frame(coeff = new_eta, se = NA),
                vcov = NA,
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG))
  } else {
    ## Calculate Cov(beta, eta) using numerical differentiation ----------------
    hN = hN_scale * N ^ ( - 1 / 2) # perturbation ------------------------------
    
    ## Calculate l(beta, eta) --------------------------------------------------
    od_loglik = mle_loglik(Y = Y,
                           offset = offset,
                           X_unval = X_unval,
                           X = X,
                           Z = Z,
                           comp_dat_val = comp_dat_val,
                           comp_dat_unval = comp_dat_unval,
                           beta = new_beta,
                           eta = new_eta, 
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
                                  Z = Z,
                                  comp_dat_val = comp_dat_val,
                                  comp_dat_unval = comp_dat_unval,
                                  beta = as.matrix(theta_pert[c(1:nrow(new_beta))]),
                                  eta = as.matrix(theta_pert[-c(1:nrow(new_beta))]), 
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
                                      Z = Z,
                                      comp_dat_val = comp_dat_val,
                                      comp_dat_unval = comp_dat_unval,
                                      beta = as.matrix(theta_pert[c(1:nrow(new_beta))]),
                                      eta = as.matrix(theta_pert[-c(1:nrow(new_beta))]), 
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
    
    ## Return final estimates and convergence information ----------------------
    return(list(coefficients = data.frame(coeff = new_beta, 
                                          se = se_beta),
                misclass_coefficients = data.frame(coeff = new_eta, 
                                                   se = se_eta),
                vcov = cov,
                converged = CONVERGED,
                se_converged = SE_CONVERGED,
                converged_msg = CONVERGED_MSG))
  }
}
