#' @title
#' Cross-validated observed-data log-likelihood
#'
#' @description
#' This function returns the value of the observed-data log-likelihood for the SMLE based on cross-validation.
#'
#' @param fold Column name with the assigned fold for cross-validation.
#' @param Y Column name with the outcome 
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = NULL} for no offset
#' @param X_unval Column name(s) with the unvalidated covariates.  If \code{X_unval} and \code{X_val} are \code{null}, all covariates are assumed to be error-free.
#' @param X_val Column name(s) with the validated covariates. If \code{X_unval} and \code{X_val} are \code{null}, all covariates are assumed to be error-free.
#' @param Z (Optional) Column name(s) with additional error-free covariates.
#' @param Validated Column name with the validation indicator. The validation indicator can be defined as \code{Validated = 1} or \code{TRUE} if the subject was validated and \code{Validated = 0} or \code{FALSE} otherwise.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{Y_unval}, \code{Y}, \code{X_unval}, \code{X_val}, \code{Z}, \code{Validated}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the covariates in the analysis model. The default assumes main effects of \code{X_val} and \code{Z} only. 
#' @param initial_lr_params Initial values for parametric model parameters. Choices include (1) \code{"Zero"} (non-informative starting values) or (2) \code{"Complete-data"} (estimated based on validated subjects only)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return scalar value of the function
#' @export

cv_smlePossum = function(fold, Y, offset = NULL, X_unval, X_val, Z = NULL, Validated = NULL, Bspline = NULL, data, theta_pred = NULL, initial_lr_params = "Zero", TOL = 1E-4, MAX_ITER = 1000) {
  if (is.null(theta_pred)) { theta_pred = c(X_val, Z) }

  num_folds = length(unique(data[, fold]))
  status = rep(TRUE, num_folds)
  msg = rep("", num_folds)
  ll = rep(NA, num_folds)
  #fold_ll = re_fold_ll = vector()
  for (i in 1:num_folds) {
    f = unique(data[, fold])[i]
    train = data[which(data[, fold] == f), ]
    train_fit = smlePossum(Y = Y, 
                           offset = offset, 
                           X_unval = X_unval, 
                           X_val = X_val, 
                           Z = Z, 
                           Validated = Validated, 
                           Bspline = Bspline, 
                           data = train, 
                           theta_pred = theta_pred, 
                           initial_lr_params = initial_lr_params, 
                           noSE = TRUE, 
                           TOL = TOL, 
                           MAX_ITER = MAX_ITER)
    status[i] = train_fit$converged
    msg[i] = train_fit$converged_msg

    if (train_fit$converged) {
      train_theta = train_fit$coeff$coeff
      train_p = train_fit$Bspline_coeff
      train_x = unique(data.frame(train[train[, Validated] == 1, X_val]))
      train_x = data.frame(train_x[order(train_x[, 1]), ])
      colnames(train_x) = X_val
      train_x = cbind(k = 1:nrow(train_x), train_x)
      train_p = merge(train_x, train_p)

      test = data[which(data[, fold] != f), ]
      test_x = unique(data.frame(test[test[, Validated] == 1, X_val]))
      test_x = data.frame(test_x[order(test_x[, 1]), ])
      colnames(test_x) = X_val
      test_x = cbind(k_ = 1:nrow(test_x), test_x)
      test_p = matrix(data = NA, nrow = nrow(test_x), ncol = length(Bspline))

      for (j in 1:nrow(test_x)) {
        x_ = test_x[j, X_val]
        bf = suppressWarnings(expr = max(which(train_x[, X_val] <= x_)))
        af = suppressWarnings(expr = min(which(train_x[, X_val] >= x_)))
        if (bf == -Inf) { bf = af }
        if (af == Inf) { af = bf }

        # X_val values immediately before/after
        x0 = train_p[bf, X_val]
        x1 = train_p[af, X_val]

        # B-spline coefficients immediately before/after
        p0 = train_p[bf, -c(1:(1 + length(X_val)))]
        p1 = train_p[af, -c(1:(1 + length(X_val)))]

        if (x1 == x0) {
          test_p[j, ] = unlist(p0)
        } else {
          test_p[j, ] = unlist((p0 * (x1 - x_) + p1 * (x_ - x0)) / (x1 - x0))
        }
      }

      # Recale columns of test_p to sum to 1
      denom = colSums(test_p)
      denom[denom == 0] = 1 # Avoid NaN error due to dividing by 0
      re_test_p = t(t(test_p) / denom)

      # Construct complete dataset  -----------------------------------
      N = nrow(test) ## total sample size (Phase I)
      n = sum(test[, Validated]) ## validation study sample size (Phase II)
      
      # Reorder so that the n validated subjects are first ------------
      test = test[order(as.numeric(test[, Validated]), decreasing = TRUE), ]
      
      # Save distinct X -------------------------------------------------
      x_obs = data.frame(unique(test[1:n, c(X_val)]))
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
      comp_dat_val = test[c(1:n), c(Y, offset, theta_pred, Bspline)]
      comp_dat_val = merge(x = comp_dat_val,
                           y = data.frame(x_obs, k = 1:m),
                           all.x = TRUE)
      comp_dat_val = comp_dat_val[, c(Y, offset, theta_pred, Bspline, "k")]
      comp_dat_val = data.matrix(comp_dat_val)
      
      # (m x n)xd vectors of each (one column per person, one row per x) --
      suppressWarnings(
        comp_dat_unval <- data.matrix(
          cbind(test[-c(1:n), c(Y, offset, setdiff(x = theta_pred, y = c(X_val)), Bspline)],
                x_obs_stacked,
                k = rep(seq(1, m), each = (N - n)))
        )
      )
      comp_dat_unval = comp_dat_unval[, c(Y, offset, theta_pred, Bspline, "k")]
      cd = rbind(comp_dat_val, comp_dat_unval)
      
      # Calculate log-likelihood -------------------------------------------
      ll_f = smle_observed_data_loglik(N = nrow(test),
                                       n = sum(test[, Validated]),
                                       Y = Y,
                                       offset = offset,
                                       X_unval = X_unval,
                                       X_val = X_val,
                                       Z = Z,
                                       Bspline = Bspline,
                                       comp_dat_all = cd, ## update
                                       theta_pred = theta_pred,
                                       theta = train_theta,
                                       p = re_test_p)
      ll[i] = ll_f
    } else {

    }
  }
  return(list(loglik = ll, status = status, msg = msg))
}

