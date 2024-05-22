#' Sieve maximum likelihood estimator (SMLE) for two-phase Poisson regression problems with covariate measurement error
#' This function returns the sieve maximum likelihood estimators (SMLE) for the Poisson regression model with covariate measurement error from Lotspeich et al. (2023+)
#'
#' @param Y Column name with the outcome 
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = NULL} for no offset
#' @param X_unval Column name(s) with the unvalidated covariates.  If \code{X_unval} and \code{X_val} are \code{null}, all covariates are assumed to be error-free.
#' @param X_val Column name(s) with the validated covariates. If \code{X_unval} and \code{X_val} are \code{null}, all covariates are assumed to be error-free.
#' @param Z (Optional) Column name(s) with additional error-free covariates.
#' @param Validated Column name with the validation indicator. The validation indicator can be defined as \code{Validated = 1} or \code{TRUE} if the subject was validated and \code{Validated = 0} or \code{FALSE} otherwise.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{Y}, \code{X_unval}, \code{X_val}, \code{Z}, \code{Validated}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the covariates in the analysis model. The default assumes main effects of \code{X_val} and \code{Z} only. 
#' @param initial_lr_params Initial values for parametric model parameters. Choices include (1) \code{"Zero"} (non-informative starting values) or (2) \code{"Complete-data"} (estimated based on validated subjects only)
#' @param h_N_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is \code{h_N_scale = 1}.
#' @param noSE Indicator for whether standard errors are desired. Defaults to \code{noSE = FALSE}.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return
#' \item{coeff}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{Bspline_coeff}{dataframe with final B-spline coefficient estimates (where applicable).}
#' \item{vcov}{variance-covariance matrix for \code{coeff} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' \item{iterations}{number of iterations completed by EM algorithm to find parameter estimates.}
#' \item{od_loglik_at_conv}{value of the observed-data log-likelihood at convergence.}
#' @export

smlePossum = function(Y, offset = NULL, X_unval, X_val, Z = NULL, Validated = NULL, Bspline = NULL, data, theta_pred = NULL, initial_lr_params = "Zero", h_N_scale = 1, noSE = FALSE, TOL = 1E-4, MAX_ITER = 1000) {
  # Prepare for algorithm -------------------------------------------
  N = nrow(data) ## total sample size (Phase I)
  n = sum(data[, Validated]) ## validation study sample size (Phase II)

  # Reorder so that the n validated subjects are first ------------
  data = data[order(as.numeric(data[, Validated]), decreasing = TRUE), ]
  # ------------------------------------------- Prepare for algorithm

  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsX = !is.null(X_unval)
  ## ------ If unvalidated variable was left blank, assume error-free
  # ----------------------------------------- Determine error setting

  # Add the B spline basis ------------------------------------------
  if (errorsX) {
    sn = ncol(data[, Bspline])
    if(0 %in% colSums(data[c(1:n), Bspline])) {
      warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)

      return(list(coeff = data.frame(coeff = NA, se = NA),
                  Bspline_coeff = NA,
                  vcov = NA,
                  converged = FALSE,
                  se_converged = FALSE,
                  converged_msg = "B-spline error",
                  iterations = 0,
                  od_loglik_at_conv = NA))
    }

  }
  # ------------------------------------------ Add the B spline basis

  
  if (is.null(theta_pred)){
    theta_pred = c(X_val, Z)
    message("Analysis model assumed main effects only.")
  }
  theta_formula = as.formula(paste0(Y, "~", paste(theta_pred, collapse = "+")))
  if (errorsX) {
    # Save distinct X -------------------------------------------------
    x_obs = data.frame(unique(data[1:n, c(X_val)]))
    x_obs = data.frame(x_obs[order(x_obs[, 1]), ])
    m = nrow(x_obs)
    x_obs_stacked = do.call(what = rbind,
                            args = replicate(n = (N - n),
                                             expr = x_obs,
                                             simplify = FALSE)
                            )
    x_obs_stacked = data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) = colnames(x_obs_stacked) = X_val

    # Save static (X*,X,Y,Z) since they don't change ---------------
    comp_dat_val = data[c(1:n), c(Y, offset, theta_pred, Bspline)]
    comp_dat_val = merge(x = comp_dat_val,
                         y = data.frame(x_obs, k = 1:m),
                         all.x = TRUE)
    comp_dat_val = comp_dat_val[, c(Y, offset, theta_pred, Bspline, "k")]
    comp_dat_val = data.matrix(comp_dat_val)

    # (m x n)xd vectors of each (one column per person, one row per x) --
    suppressWarnings(
      comp_dat_unval <- data.matrix(
        cbind(data[-c(1:n), c(Y, offset, setdiff(x = theta_pred, y = c(X_val)), Bspline)],
              x_obs_stacked,
              k = rep(seq(1, m), each = (N - n)))
        )
      )
    comp_dat_unval = comp_dat_unval[, c(Y, offset, theta_pred, Bspline, "k")]

    comp_dat_all = rbind(comp_dat_val, comp_dat_unval)

    # Initialize B-spline coefficients {p_kj}  ------------
    ## Numerators sum B(Xi*) over k = 1,...,m -------------
    ## Save as p_val_num for updates ----------------------
    ## (contributions don't change) -----------------------
    p_val_num = rowsum(x = comp_dat_val[, Bspline],
                       group = comp_dat_val[, "k"],
                       reorder = TRUE)
    prev_p = p0 =  t(t(p_val_num) / colSums(p_val_num))
  }
  theta_design_mat = cbind(int = 1, comp_dat_all[, theta_pred])

  # Initialize parameter values -------------------------------------
  ## theta, gamma ---------------------------------------------------
  if(!(initial_lr_params %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    initial_lr_params = "Zero"
  }

  if(initial_lr_params == "Zero") {
    prev_theta = theta0 = matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
  } else if(initial_lr_params == "Complete-data") {
    if (!is.null(offset)) {
      prev_theta = theta0 = matrix(glm(formula = theta_formula,
                                       family = "poisson",
                                       data = data.frame(data[c(1:n), ]),
                                       offset = log(data[c(1:n), offset]))$coefficients,
                                   ncol = 1)
    } else {
      prev_theta = theta0 = matrix(glm(formula = theta_formula,
                                       family = "poisson",
                                       data = data.frame(data[c(1:n), ]))$coefficients,
                                   ncol = 1)
    }
  }

  # Set parameters for algorithm convergence --------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  
  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta = as.numeric((theta_design_mat[-c(1:n), ] %*% prev_theta))
    lambda = exp(mu_theta)
    if (!is.null(offset)) {
       lambda = comp_dat_unval[, offset] * lambda
    }
    pY_X = dpois(x = comp_dat_unval[, Y], 
                 lambda = lambda)
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(X|X*) -------------------------------------------------------
    if (errorsX) {
      ### p_kj ----------------------------------------------------------
      ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
      ### multiply by the B-spline terms
      pX = prev_p[rep(seq(1, m), each = (N - n)), ] * comp_dat_unval[, Bspline]
      ### ---------------------------------------------------------- p_kj
    }
    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------
    if (errorsX) {
      ### P(Y|X,Z)p_kjB(X*) -------------------------------------------
      psi_num = c(pY_X) * pX
      ### Update denominator ------------------------------------------
      #### Sum up all rows per id (e.g. sum over xk) ------------------
      psi_denom = rowsum(x = psi_num,
                         group = rep(seq(1, (N - n)), times = m))
      #### Then sum over the sn splines -------------------------------
      psi_denom = rowSums(psi_denom)
      #### Avoid NaN resulting from dividing by 0 ---------------------
      psi_denom[psi_denom == 0] = 1
      ### And divide them! --------------------------------------------
      psi_t = psi_num / psi_denom
      ### Update the w_kyi for unvalidated subjects -------------------
      ### by summing across the splines/ columns of psi_t -------------
      w_t = rowSums(psi_t)
    }
    ### ----------------------------- Estimate conditional expectations
    # ---------------------------------------------------------- E Step
    ###################################################################

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta using weighted Poisson regression -----------------
    ### Gradient ------------------------------------------------------
    w_t = c(rep(1, n), w_t) #w_t = lengthenWT(w_t, n)
    # muVector = calculateMu(theta_design_mat, prev_theta)
    # gradient_theta = calculateGradient(w_t, n, theta_design_mat, comp_dat_all[, Y], muVector)
    # hessian_theta = calculateHessian(theta_design_mat, w_t, muVector, n, mus_theta);
    # ### ------------------------------------------------------ Gradient
    # ### Hessian -------------------------------------------------------
    # new_theta = tryCatch(expr = prev_theta - (solve(hessian_theta) %*% gradient_theta),
    #                      error = function(err) { matrix(NA, nrow = nrow(prev_theta)) })
    # if (any(is.na(new_theta)))
    # {
    if (!is.null(offset)) {
      suppressWarnings(
        new_theta <- matrix(data = glm(formula = theta_formula,
                                       family = "poisson",
                                       data = data.frame(comp_dat_all),
                                       weights = w_t, 
                                       offset = log(comp_dat_all[, offset]))$coefficients,
                            ncol = 1)
      )
    } else {
      suppressWarnings(
        new_theta <- matrix(data = glm(formula = theta_formula,
                                       family = "poisson",
                                       data = data.frame(comp_dat_all),
                                       weights = w_t)$coefficients,
                            ncol = 1)
      )
    }
    #}

    ### Check for convergence -----------------------------------------
    theta_conv = abs(new_theta - prev_theta) < TOL

    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    if (errorsX) {
      ### Update numerators by summing u_t over i = 1, ..., N ---------
      new_p_num = p_val_num +
        rowsum(x = psi_t,
               group = rep(x = seq(1, m),
                           each = (N - n)),
               reorder = TRUE)
      new_p = t(t(new_p_num) / colSums(new_p_num))
      ### Check for convergence ---------------------------------------
      p_conv <- abs(new_p - prev_p) < TOL
    } else {
      p_conv <- TRUE
    }
    ## -------------------------------------------------- Update {p_kj}
    # ---------------------------------------------------------- M Step
    ###################################################################
    # Check for convergence -------------------------------------------
    all_conv = c(theta_conv, p_conv)
    CONVERGED = mean(all_conv) == 1

    # Update values for next iteration  -------------------------------
    it = it + 1
    prev_theta = new_theta
    if (errorsX) {
      prev_p = new_p
    }
    #  ------------------------------- Update values for next iteration
  }
  rownames(new_theta) <- c("Intercept", theta_pred)

  if(!CONVERGED) {
    if(it > MAX_ITER) {
      CONVERGED_MSG = "MAX_ITER reached"
    }

    return(list(coeff = data.frame(coeff = NA, se = NA),
                outcome_err_coeff = data.frame(coeff = NA, se = NA),
                Bspline_coeff = NA,
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

  # ---------------------------------------------- Estimate theta using EM
  if(noSE){
    if (!errorsX) {
      new_p = p_val_num = matrix(data = NA, 
                                 nrow = 1,
                                 ncol = 1)
    }

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta = smle_observed_data_loglik(N = N,
                                                n = n,
                                                Y = Y,
                                                offset = offset,
                                                X_unval = X_unval,
                                                X_val = X_val,
                                                Z = Z,
                                                Bspline = Bspline,
                                                comp_dat_all = comp_dat_all,
                                                theta_pred = theta_pred,
                                                theta = new_theta,
                                                p = new_p)

    return(list(coeff = data.frame(coeff = new_theta, se = NA),
                Bspline_coeff = cbind(k = 1:m, new_p),
                vcov = NA,
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_theta))
  } else {
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_N = h_N_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    if (!errorsX) {
      new_p = p_val_num = matrix(NA, 
                                 nrow = 1, 
                                 ncol = 1)
    }

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta = smle_observed_data_loglik(N = N,
      n = n,
      Y = Y,
      X_unval = X_unval,
      X_val = X_val,
      Z = Z,
      Bspline = Bspline,
      comp_dat_all = comp_dat_all,
      theta_pred = theta_pred,
      theta = new_theta,
      p = new_p)

    I_theta <- matrix(od_loglik_theta, 
                      nrow = nrow(new_theta), 
                      ncol = nrow(new_theta))

    single_pert_theta <- sapply(X = seq(1, ncol(I_theta)),
      FUN = pl_theta,
      theta = new_theta,
      h_N = h_N,
      n = n,
      N = N,
      Y = Y,
      offset = offset, 
      X_unval = X_unval,
      X_val = X_val,
      Z = Z,
      Bspline = Bspline,
      comp_dat_all = comp_dat_all,
      theta_pred = theta_pred,
      p0 = new_p,
      p_val_num = p_val_num,
      TOL = TOL,
      MAX_ITER = MAX_ITER)

    if (any(is.na(single_pert_theta))) {
      I_theta <- matrix(NA, nrow = nrow(new_theta), ncol = nrow(new_theta))
      SE_CONVERGED <- FALSE
    } else {
      spt_wide <- matrix(rep(c(single_pert_theta), times = ncol(I_theta)),
       ncol = ncol(I_theta),
       byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
      I_theta <- I_theta - spt_wide - t(spt_wide)
      SE_CONVERGED <- TRUE
    }

    for (c in 1:ncol(I_theta)) {
      pert_theta <- new_theta
      pert_theta[c] <- pert_theta[c] + h_N
      double_pert_theta <- sapply(X = seq(c, ncol(I_theta)),
        FUN = pl_theta,
        theta = pert_theta,
        h_N = h_N,
        n = n,
        N = N,
        Y = Y,
        offset = offset,
        X_unval = X_unval,
        X_val = X_val,
        Z = Z,
        Bspline = Bspline,
        comp_dat_all = comp_dat_all,
        theta_pred = theta_pred,
        p0 = new_p,
        p_val_num = p_val_num,
        MAX_ITER = MAX_ITER,
        TOL = TOL)

      dpt <- matrix(0, nrow = nrow(I_theta), ncol = ncol(I_theta))
      dpt[c,c] <- double_pert_theta[1] #Put double on the diagonal
      if(c < ncol(I_theta)) {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] <- dpt[-(1:c), c] <- double_pert_theta[-1]
      }

      I_theta <- I_theta + dpt
    }

    I_theta <- h_N ^ (- 2) * I_theta

    cov_theta <- tryCatch(expr = - solve(I_theta),
      error = function(err) {
        matrix(NA, nrow = nrow(I_theta), ncol = ncol(I_theta)) }
      )
    # ------------------------- Estimate Cov(theta) using profile likelihood

    se_theta <- tryCatch(expr = sqrt(diag(cov_theta)),
      warning = function(w) {
        matrix(NA, nrow = nrow(prev_theta))}
      )
    SE_CONVERGED <- any(is.na(se_theta))

    return(list(coeff = data.frame(coeff = new_theta, se = se_theta),
                Bspline_coeff = cbind(X = x_obs, new_p),
                vcov = cov_theta,
                converged = CONVERGED,
                se_converged = SE_CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_theta))
  }
}

