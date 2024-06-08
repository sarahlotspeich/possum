#' Sieve maximum likelihood estimator (SMLE) for two-phase Poisson regression problems with covariate measurement error
#' This function returns the sieve maximum likelihood estimators (SMLE) for the Poisson regression model with covariate measurement error from Lotspeich et al. (2023+)
#'
#' @param Y Column name with the outcome 
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = NULL} for no offset.
#' @param X_unval Column name(s) with the unvalidated covariates.  
#' @param X Column name(s) with the validated covariates.
#' @param Z (Optional) Column name(s) with additional error-free covariates.
#' @param data A dataframe with one row per subject containing columns: \code{Y}, \code{X_unval}, \code{X}, \code{Z}, \code{Validated}, and \code{Bspline}.
#' @param beta_pred Vector of columns in \code{data} that pertain to the covariates in the analysis model. The default assumes main effects of \code{X} and \code{Z} only. 
#' @param beta0 Initial values for parametric model parameters. Choices include (1) \code{"Zero"} (non-informative starting values, the default) or (2) \code{"Complete-data"} (estimated based on validated subjects only)
#' @param h_N_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is \code{h_N_scale = 1}.
#' @param noSE Indicator for whether standard errors are desired. Defaults to \code{noSE = FALSE}.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return
#' \item{coeff}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{ppv}{final positive predictive value estimate.}
#' \item{vcov}{variance-covariance matrix for \code{coeff} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' \item{iterations}{number of iterations completed by EM algorithm to find parameter estimates.}
#' \item{od_loglik_at_conv}{value of the observed-data log-likelihood at convergence.}
#' @export

smlePossum = function(Y, offset = NULL, X_unval, X, Z = NULL, data, beta_pred = NULL, beta0 = "Zero", h_N_scale = 1, noSE = FALSE, TOL = 1E-4, MAX_ITER = 1000) {
  # Prepare for algorithm ------------------------------------------------------
  data[, "Validated"] = as.numeric(!is.na(data[, X])) ## validation indicator
  N = nrow(data) ## total sample size (Phase I)
  n = sum(data[, "Validated"]) ## validation study sample size (Phase II)

  ## Reorder so that the n validated subjects are first ------------------------
  data = data[order(as.numeric(data[, "Validated"]), decreasing = TRUE), ]
  
  ## Create row numbers --------------------------------------------------------
  data[, "row_num"] = 1:N
  ##############################################################################
  ## Define formula for analysis model -----------------------------------------
  ### If beta_pred not specified, default to main effects only -----------------
  if (is.null(beta_pred)){
    beta_pred = c(X, Z)
  }
  beta_formula = as.formula(paste0(Y, "~", paste(beta_pred, collapse = "+")))
  ##############################################################################
  ## Save static (X*,X,Y,Z) for validated rows since they don't change ---------
  comp_dat_val = data.matrix(data[c(1:n), c(Y, offset, X_unval, beta_pred, "row_num")])
  
  ## Create augmented (X*,x,Y,Z) for unvalidated rows --------------------------
  ### First (N-n) rows assume X = 0 --------------------------------------------
  comp_dat0 = data.matrix(data[-c(1:n), c(Y, offset, X_unval, beta_pred, "row_num")])
  comp_dat0[, X] = 0 
  ### Last (N-n) rows assume X = 1 ---------------------------------------------
  comp_dat1 = data.matrix(data[-c(1:n), c(Y, offset, X_unval, beta_pred, "row_num")])
  comp_dat1[, X] = 1
  ### Put them together --------------------------------------------------------
  comp_dat_unval = rbind(comp_dat0, 
                         comp_dat1)
  
  ## Create augmented "complete" dataset of validated and unvalidated ----------
  comp_dat_all = rbind(comp_dat_val, comp_dat_unval)
  ##############################################################################
  ## Initialize analysis model parameters (beta) -------------------------------
  ### If invalid beta0 specified, default to beta0 = "Zero" --------------------
  if(!(beta0 %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    beta0 = "Zero"
  }
  ### Set initial values for beta ----------------------------------------------
  if(beta0 == "Zero") {
    prev_beta = beta0 = matrix(data = 0, 
                               nrow = (length(beta_pred) + 1), 
                               ncol = 1)
  } else if(beta0 == "Complete-data") {
    if (!is.null(offset)) {
      prev_beta = beta0 = matrix(glm(formula = beta_formula,
                                     family = "poisson",
                                     data = data.frame(data[c(1:n), ]),
                                     offset = log(data[c(1:n), offset]))$coefficients,
                                 ncol = 1)
    } else {
      prev_beta = beta0 = matrix(glm(formula = beta_formula,
                                     family = "poisson",
                                     data = data.frame(data[c(1:n), ]))$coefficients,
                                 ncol = 1)
    }
  }
  ##############################################################################
  ## Initialize positive predictive value (PPV) --------------------------------
  TP = sum(comp_dat_val[, X] == 1 & comp_dat_val[, X_unval] == 1) ### number of true positives in validated subset
  FP = sum(comp_dat_val[, X] == 0 & comp_dat_val[, X_unval] == 1)  ### number of false positives in validated subset
  prev_ppv = ppv0 = TP / (TP + FP)
  # ------------------------------------------------------ Prepare for algorithm
  ##############################################################################
  # Estimate beta and ppv using EM algorithm -----------------------------------
  ## Set parameters for algorithm convergence ----------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  ## Initialize components of ppv update ---------------------------------------
  ### Numerator sums over Xi* among validated rows with Xi = 1 -----------------
  ppv_num0 = sum(comp_dat_val[comp_dat_val[, X] == 1, X_unval])
  ### Denominator sums over Xi* among all validated and unvalidated rows -------
  ppv_denom0 = sum(data[, X_unval])
  ## Begin algorithm -----------------------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step -------------------------------------------------------------------
    ## Update the psi_xi = P(X=x|Yi,Xi*,Z) for unvalidated subjects ------------
    ### Analysis model: P(Y|X) -------------------------------------------------
    #### mu = beta0 + beta1X + beta2Z + ... 
    mu_beta = as.numeric(cbind(int = 1, comp_dat_unval[, beta_pred]) %*% prev_beta)
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
    ### Misclassification mechanism: P(X|X*) -----------------------------------
    #### P(X|X*=1) = ppv ^ X * (1 - X) -----------------------------------------
    pXgivXstar = prev_ppv ^ comp_dat_unval[, X] * (1 - prev_ppv) ^ (1 - comp_dat_unval[, X])
    #### P(X|X*=0) = I(X=0) since false negatives aren't possible --------------
    pXgivXstar[comp_dat_unval[, X_unval] == 0] = as.numeric(comp_dat_unval[comp_dat_unval[, X_unval] == 0, X] == 0)
    ############################################################################
    ## Estimate conditional expectations --------------------------------------
    ### Update numerator -------------------------------------------------------
    #### P(Y|X,Z)P(X|X*) -------------------------------------------------------
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
    if (!is.null(offset)) {
      new_beta = suppressWarnings(
        matrix(data = glm(formula = beta_formula,
                          family = "poisson",
                          data = data.frame(comp_dat_all),
                          weights = psi_aug, 
                          offset = log(comp_dat_all[, offset]))$coefficients,
               ncol = 1)
      )
    } else {
      new_beta = suppressWarnings(
        matrix(data = glm(formula = beta_formula,
                          family = "poisson",
                          data = data.frame(comp_dat_all),
                          weights = psi_aug)$coefficients,
               ncol = 1)
      )
    }
    ## Check for beta convergence ----------------------------------------------
    beta_conv = abs(new_beta - prev_beta) < TOL
    ############################################################################
    ## Update ppv --------------------------------------------------------------
    ### Update numerator by adding sum over X* x psi_1i for unvalidated --------
    ppv_num = ppv_num0 + 
      sum(comp_dat_unval[-c(1:(N-n)), X_unval] * psi[-c(1:(N-n))])
    ### Divide them ------------------------------------------------------------
    new_ppv = ppv_num / ppv_denom0
    ### Check for ppv convergence ----------------------------------------------
    ppv_conv = abs(new_ppv - prev_ppv) < TOL
    ############################################################################
    # Check for global convergence ---------------------------------------------
    all_conv = c(beta_conv, ppv_conv)
    CONVERGED = mean(all_conv) == 1
    # Update values for next iteration  ----------------------------------------
    it = it + 1
    prev_beta = new_beta
    prev_ppv = new_ppv
  }
  rownames(new_beta) <- c("Intercept", beta_pred)

  if(!CONVERGED) {
    if(it > MAX_ITER) {
      CONVERGED_MSG = "MAX_ITER reached"
    }

    return(list(coeff = data.frame(coeff = NA, se = NA),
                ppv = NA,
                vcov = NA,
                converged = FALSE,
                se_converged = NA,
                converged_msg = "MAX_ITER reached",
                iterations = it,
                od_loglik_at_conv = NA))
  }

  if(CONVERGED) { 
    CONVERGED_MSG = "Converged" 
  }

  # ---------------------------------------------- Estimate beta using EM
  if(noSE){
    ## Calculate l(beta) --------------------------------------------------
    od_loglik_beta = smle_loglik(Y = Y,
                                 offset = offset,
                                 X_unval = X_unval,
                                 X = X,
                                 Z = Z,
                                 comp_dat_val = comp_dat_val,
                                 comp_dat_unval = comp_dat_unval,
                                 beta_pred = beta_pred,
                                 beta = new_beta,
                                 ppv = new_ppv)

    return(list(coeff = data.frame(coeff = new_beta, se = NA),
                ppv = new_ppv,
                vcov = NA,
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_beta))
  } else {
    # Estimate Cov(beta) using profile likelihood ------------------------------
    h_N = h_N_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    ## Calculate pl(beta) ------------------------------------------------------
    od_loglik_beta = smle_loglik(Y = Y,
                                  offset = offset,                                          
                                  X_unval = X_unval,
                                  X = X,
                                  Z = Z,
                                  comp_dat_val = comp_dat_val,
                                  comp_dat_unval = comp_dat_unval,
                                  beta_pred = beta_pred,
                                  beta = new_beta,
                                  ppv = new_ppv)
    
    ## Calculate I(beta) using profile likelihood method -----------------------
    I_beta = matrix(data = od_loglik_beta, 
                    nrow = nrow(new_beta), 
                    ncol = nrow(new_beta))
  
    ### Profile log-likelihoods after single perturbations of beta -------------
    single_pert_beta = vector(length = ncol(I_beta))
    for (i in 1:length(single_pert_beta)) {
      single_pert_beta[i] = pl_beta(k = i, 
                                    beta = new_beta, 
                                    h_N = h_N, 
                                    Y = Y, 
                                    offset = offset, 
                                    X_unval = X_unval, 
                                    X = X, 
                                    Z = Z, 
                                    comp_dat_val = comp_dat_val, 
                                    comp_dat_unval = comp_dat_unval, 
                                    beta_pred = beta_pred, 
                                    ppv0 = ppv0, 
                                    ppv_num0 = ppv_num0, 
                                    ppv_denom0 = ppv_denom0, 
                                    TOL = TOL, 
                                    MAX_ITER = MAX_ITER)
    }
    
    #### Check for any elements that didn't converge ---------------------------
    if (any(is.na(single_pert_beta))) {
      I_beta = matrix(data = NA, 
                      nrow = nrow(new_beta), 
                      ncol = nrow(new_beta))
    } else {
      ### Create wide version of single perturbations --------------------------
      spt_wide = matrix(data = rep(c(single_pert_beta), 
                                   times = ncol(I_beta)),
                        ncol = ncol(I_beta),
                        byrow = FALSE)
      #for the each kth row of single_pert_beta add to the kth row / kth column of I_beta
      I_beta = I_beta - spt_wide - t(spt_wide)
      
      ### Profile log-likelihoods after double perturbations of beta -----------
      #### Create matrix of zeros for them (to be added to I_beta) -------------
      dpb = matrix(data = 0, 
                   nrow = nrow(I_beta), 
                   ncol = ncol(I_beta))
      #### Loop over the rows and columns to fill dpb --------------------------
      for (r in 1:nrow(dpb)) {
        #### First perturb the rth element in beta 
        pert_beta = new_beta
        pert_beta[r] = pert_beta[r] + h_N
        #### Then loop over further perturbing all elements in beta
        for (c in r:ncol(dpb)) {
          dpb[r, c] = dpb[c, r] = pl_beta(k = c, 
                                          beta = pert_beta, 
                                          h_N = h_N, 
                                          Y = Y, 
                                          offset = offset, 
                                          X_unval = X_unval, 
                                          X = X, 
                                          Z = Z, 
                                          comp_dat_val = comp_dat_val, 
                                          comp_dat_unval = comp_dat_unval, 
                                          beta_pred = beta_pred, 
                                          ppv0 = ppv0, 
                                          ppv_num0 = ppv_num0, 
                                          ppv_denom0 = ppv_denom0, 
                                          TOL = TOL, 
                                          MAX_ITER = MAX_ITER)
          
        }
      }
      
      ### Add double perturbations to matrix of singles ------------------------
      I_beta = I_beta + dpb
    }

    I_beta = - h_N ^ (- 2) * I_beta

    cov_beta = tryCatch(expr = solve(I_beta),
      error = function(err) {
        matrix(NA, nrow = nrow(I_beta), ncol = ncol(I_beta)) }
      )
    # ------------------------- Estimate Cov(beta) using profile likelihood

    se_beta <- tryCatch(expr = sqrt(diag(cov_beta)),
      warning = function(w) {
        matrix(NA, nrow = nrow(prev_beta))}
      )
    SE_CONVERGED <- any(is.na(se_beta))

    return(list(coeff = data.frame(coeff = new_beta, se = se_beta),
                ppv = new_ppv,
                vcov = cov_beta,
                converged = CONVERGED,
                se_converged = SE_CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_beta))
  }
}

