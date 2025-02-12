#' Observed-data log-likelihood for the sieve maximum likelihood estimator (SMLE)
#'
#' This function returns the value of the observed-data log-likelihood (equation (#) in Lotspeich et al. (2023+))
#' for a given dataset and parameter values `beta` and `p`.
#
#'
#' @param Y Column name with the outcome
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated covariates
#' @param X Column(s) with the validated covariates
#' @param Z (Optional) Column(s) with additional error-free covariates
#' @param comp_dat_val Dataset containing rows for validated subjects' data (a matrix)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param beta Parameters for the analysis model (a column vector)
#' @param eta Parameters for the misclassification model (a column vector)
#' @param noFN logical, if \code{noFN = FALSE} (the default), then it is assumed that there can be both false positives and false negatives in the error-prone exposure. If \code{noFN = TRUE}, the error mechanism is restricted to only false positives.
#' @return Scalar value of the function

mle_loglik = function(Y = NULL, offset = NULL, X_unval = NULL, X = NULL, Z = NULL, comp_dat_val, comp_dat_unval, beta, eta, noFN) {
  ##############################################################################
  # Save useful constants ------------------------------------------------------
  N = max(comp_dat_unval[, "row_num"]) ## total sample size
  n = max(comp_dat_val[, "row_num"]) ## validation sub-sample size

  # Create combined dataset of validated and unvalidated -----------------------
  comp_dat_all = rbind(comp_dat_val,
                       comp_dat_unval)
  ##############################################################################
  # Calculate probabilities ----------------------------------------------------
  ## Analysis model: P(Y|X) ----------------------------------------------------
  ### mu = beta0 + beta1X + beta2Z + ...
  mu_beta = as.numeric(cbind(int = 1, comp_dat_all[, c(X, Z)]) %*% beta)
  ### lambda = exp(beta0 + beta1X + beta2Z + ... )
  lambda = exp(mu_beta)
  ### If offset specified, lambda = offset x exp(beta0 + beta1X + beta2Z + ... )
  if (!is.null(offset)) {
    lambda = comp_dat_all[, offset] * lambda
  }
  ### Calculate P(Y|X) from Poisson distribution
  pYgivX = dpois(x = comp_dat_all[, Y],
                 lambda = lambda)
  ##############################################################################
  ## Misclassification mechanism: P(X|X*) --------------------------------------
  if (noFN) { #### If one-sided errors, logistic regression on just X*=1 -----
    #### mu = eta0 + eta1Z + ...
    mu_eta = as.numeric(cbind(int = 1, comp_dat_all[, Z]) %*% eta)
    #### Calculate P(X|X*=1,Z) from Bernoulli distribution -------------------
    pXgivXstar = dbinom(x = comp_dat_all[, X],
                        size = 1,
                        prob = 1 / (1 + exp(- mu_eta)))
    #### Force P(X=0|X*=0,Z)=1 and P(X=1|X*=0,Z)=0 for all Z -----------------
    pXgivXstar[which(comp_dat_all[, X_unval] == 0 & comp_dat_all[, X] == 0)] = 1
    pXgivXstar[which(comp_dat_all[, X_unval] == 0 & comp_dat_all[, X] == 1)] = 0
  } else { #### If two-sided errors, logistic regression on all rows ---------
    #### mu = eta0 + eta1X* + eta2Z + ...
    mu_eta = as.numeric(cbind(int = 1, comp_dat_all[, c(X_unval, Z)]) %*% eta)
    #### Calculate P(X|X*,Z) from Bernoulli distribution ---------------------
    pXgivXstar = dbinom(x = comp_dat_all[, X],
                        size = 1,
                        prob = 1 / (1 + exp(- mu_eta)))
  }
  ##############################################################################
  ## Joint conditional: P(Y,X|X*) ----------------------------------------------
  pYXgivXstar = pYgivX * pXgivXstar
  ##############################################################################
  # Log-likelihood contribution of validated observations  ---------------------
  ## Sum over log{P(Y|X)} for validated observations (first n rows)
  log_pYgivX = log(pYgivX[1:n])
  log_pYgivX[log_pYgivX == -Inf] = 0
  return_loglik = sum(log_pYgivX)
  ## Sum over log{P(X|X*)}
  log_pXgivXstar = log(pXgivXstar[1:n])
  log_pXgivXstar[log_pXgivXstar == -Inf] = 0
  return_loglik = return_loglik + sum(log_pXgivXstar)
  ##############################################################################
  # Log-likelihood contribution of unvalidated observations  -------------------
  ## Marginalize out X: P(Y|X*) ------------------------------------------------
  pYXgivXstar_unval = pYXgivXstar[-c(1:n)] ### remove validated observations (first n rows)
  pYgivXstar_unval = pYXgivXstar_unval[1:(N - n)] + ### sum over P(Y,X=0|X*) (first N - n rows)
    pYXgivXstar_unval[-c(1:(N - n))] ### and P(Y,X=1|X*) (last N - n rows)
  ## Sum over log{P(Y|X*)}
  log_pYgivXstar_unval = log(pYgivXstar_unval)
  log_pYgivXstar_unval[log_pYgivXstar_unval == -Inf] = 0
  return_loglik = return_loglik + sum(log_pYgivXstar_unval)
  return(return_loglik)
}

mle_loglik_nd = function(beta_eta, dim_beta,
                         Y = NULL, offset = NULL, X_unval = NULL, X = NULL, Z = NULL,
                         comp_dat_val, comp_dat_unval, noFN) {
  ##############################################################################
  # Save useful constants ------------------------------------------------------

  ## split parameters into analysis and error model parameters
  beta <- beta_eta[1:dim_beta]
  eta <- beta_eta[-c(1:dim_beta)]

  ## save sample sizes
  N = max(comp_dat_unval[, "row_num"]) ## total sample size
  n = max(comp_dat_val[, "row_num"]) ## validation sub-sample size

  # Create combined dataset of validated and unvalidated -----------------------
  comp_dat_all = rbind(comp_dat_val,
                       comp_dat_unval)
  ##############################################################################
  # Calculate probabilities ----------------------------------------------------
  ## Analysis model: P(Y|X) ----------------------------------------------------
  ### mu = beta0 + beta1X + beta2Z + ...
  mu_beta = as.numeric(cbind(int = 1, comp_dat_all[, c(X, Z)]) %*% beta)
  ### lambda = exp(beta0 + beta1X + beta2Z + ... )
  lambda = exp(mu_beta)
  ### If offset specified, lambda = offset x exp(beta0 + beta1X + beta2Z + ... )
  if (!is.null(offset)) {
    lambda = comp_dat_all[, offset] * lambda
  }
  ### Calculate P(Y|X) from Poisson distribution
  pYgivX = dpois(x = comp_dat_all[, Y],
                 lambda = lambda)
  ##############################################################################
  ## Misclassification mechanism: P(X|X*) --------------------------------------
  if (noFN) { #### If one-sided errors, logistic regression on just X*=1 -----
    #### mu = eta0 + eta1Z + ...
    mu_eta = as.numeric(cbind(int = 1, comp_dat_all[, Z]) %*% eta)
    #### Calculate P(X|X*=1,Z) from Bernoulli distribution -------------------
    pXgivXstar = dbinom(x = comp_dat_all[, X],
                        size = 1,
                        prob = 1 / (1 + exp(- mu_eta)))
    #### Force P(X=0|X*=0,Z)=1 and P(X=1|X*=0,Z)=0 for all Z -----------------
    pXgivXstar[which(comp_dat_all[, X_unval] == 0 & comp_dat_all[, X] == 0)] = 1
    pXgivXstar[which(comp_dat_all[, X_unval] == 0 & comp_dat_all[, X] == 1)] = 0
  } else { #### If two-sided errors, logistic regression on all rows ---------
    #### mu = eta0 + eta1X* + eta2Z + ...
    mu_eta = as.numeric(cbind(int = 1, comp_dat_all[, c(X_unval, Z)]) %*% eta)
    #### Calculate P(X|X*,Z) from Bernoulli distribution ---------------------
    pXgivXstar = dbinom(x = comp_dat_all[, X],
                        size = 1,
                        prob = 1 / (1 + exp(- mu_eta)))
  }
  ##############################################################################
  ## Joint conditional: P(Y,X|X*) ----------------------------------------------
  pYXgivXstar = pYgivX * pXgivXstar
  ##############################################################################
  # Log-likelihood contribution of validated observations  ---------------------
  ## Sum over log{P(Y|X)} for validated observations (first n rows)
  log_pYgivX = log(pYgivX[1:n])
  log_pYgivX[log_pYgivX == -Inf] = 0
  return_loglik = sum(log_pYgivX)
  ## Sum over log{P(X|X*)}
  log_pXgivXstar = log(pXgivXstar[1:n])
  log_pXgivXstar[log_pXgivXstar == -Inf] = 0
  return_loglik = return_loglik + sum(log_pXgivXstar)
  ##############################################################################
  # Log-likelihood contribution of unvalidated observations  -------------------
  ## Marginalize out X: P(Y|X*) ------------------------------------------------
  pYXgivXstar_unval = pYXgivXstar[-c(1:n)] ### remove validated observations (first n rows)
  pYgivXstar_unval = pYXgivXstar_unval[1:(N - n)] + ### sum over P(Y,X=0|X*) (first N - n rows)
    pYXgivXstar_unval[-c(1:(N - n))] ### and P(Y,X=1|X*) (last N - n rows)
  ## Sum over log{P(Y|X*)}
  log_pYgivXstar_unval = log(pYgivXstar_unval)
  log_pYgivXstar_unval[log_pYgivXstar_unval == -Inf] = 0
  return_loglik = return_loglik + sum(log_pYgivXstar_unval)
  return(return_loglik)
}

loglik_mat = function(beta_eta,
                      Y_name, X_name,
                      Z_name = NULL, Xstar_name,
                      Q_name,
                      offset_name = NULL,
                      data,
                      noFN = FALSE,
                      verbose = FALSE) {
  # Save useful constants
  N = nrow(data) ## Phase I sample size
  n = sum(data[, Q_name]) ## Phase II sample size

  # Reorder data to put queried rows first
  data = data[order(data[, Q_name], decreasing = TRUE), ]

  # Create matrix of complete data
  if(!is.null(Z_name)){ #case with covariates
    if (n < N) {
      queried_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Z_name, Xstar_name, offset_name)])
      unqueried_data = rbind(
        cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 0, data[-c(1:n), c(Z_name, Xstar_name, offset_name)]),
        cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 1, data[-c(1:n), c(Z_name, Xstar_name, offset_name)])
      )
      colnames(unqueried_data) = c("id", Y_name, X_name, Z_name, Xstar_name, offset_name)
      if (noFN) {
        ## If false negatives aren't possible, delete rows from the unqueried complete data
        ### where X = 1 and X* = 0 (false negative)
        unqueried_data = unqueried_data[!(unqueried_data[, X_name] == 1 & unqueried_data[, Xstar_name] == 0), ]
      }
      complete_data = data.matrix(rbind(queried_data, unqueried_data))
    } else {
      complete_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Z_name, Xstar_name, offset_name)])
    }
  } else{ #case without covariates
    if (n < N) {
      queried_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Z_name, Xstar_name, offset_name)])
      unqueried_data = rbind(
        cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 0, data[-c(1:n), c(Xstar_name, offset_name)]),
        cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 1, data[-c(1:n), c(Xstar_name, offset_name)])
      )
      colnames(unqueried_data) = c("id", Y_name, X_name, Xstar_name, offset_name)
      if (noFN) {
        ## If false negatives aren't possible, delete rows from the unqueried complete data
        ### where X = 1 and X* = 0 (false negative)
        unqueried_data = unqueried_data[!(unqueried_data[, X_name] == 1 & unqueried_data[, Xstar_name] == 0), ]
      }
      complete_data = data.matrix(rbind(queried_data, unqueried_data))
    } else {
      complete_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Xstar_name, offset_name)])
      }
  }

  # Compute log-likelihood
  if(!is.null(Z_name)){ ## P(Y|X,Z) from Poisson distribution
    lambdaY = exp(beta_eta[1] + beta_eta[2] * complete_data[, X_name] + beta_eta[3] * complete_data[, Z_name])
    if(!is.null(offset_name)){ #has offset_name
      lambdaY = complete_data[, offset_name] * lambdaY
    }
  } else{ ## P(Y|X) from Poisson distribution
    lambdaY = exp(beta_eta[1] + beta_eta[2] * complete_data[, X_name])
    if(!is.null(offset_name)){ #has offset
      lambdaY = complete_data[, offset_name] * lambdaY
    }
  }
  ### Dazzle fix: replace y with data[, Y_name]
  pYgivXZ = dpois(x = complete_data[, Y_name], lambda = lambdaY)

  ## P(X|X*,Z) from Bernoulli distribution
  ### or mixture if noFN = TRUE
  if (noFN) {
    if(!is.null(Z_name)){ ## P(X|X*,Z) from Bernoulli distribution
      expit_XgivXstarZ = 1 / (1 + exp(-(beta_eta[4] + beta_eta[5] * complete_data[, Z_name])))
    } else{ ## P(X|X*) from Bernoulli distribution
      expit_XgivXstarZ = 1 / (1 + exp(-(beta_eta[3])))
    }

    ## P(X|X*,Z) from Bernoulli distribution
    pXgivXstarZ = expit_XgivXstarZ ^ complete_data[, X_name] * (1 - expit_XgivXstarZ) ^ (1 - complete_data[, X_name])

    ## But if X* = 0, replace with point mass
    # pXgivXstarZ[complete_data[, Xstar_name] == 0 & complete_data[, X_name] == 0] = 1 ### P(X=0|X*=0) = 1
    # pXgivXstarZ[complete_data[, Xstar_name] == 0 & complete_data[, X_name] == 1] = 0 ### P(X=1|X*=0) = 0
  } else {
    if(!is.null(Z_name)){ ## P(X|X*,Z) from Bernoulli distribution
      expit_XgivXstarZ = 1 / (1 + exp(-(beta_eta[4] + beta_eta[5] * complete_data[, Xstar_name] + beta_eta[6] * complete_data[, Z_name])))
    } else{ ## P(X|X*) from Bernoulli distribution
      expit_XgivXstarZ = 1 / (1 + exp(-(beta_eta[3] + beta_eta[4] * complete_data[, Xstar_name])))
    }

    ## P(X|X*,Z) from Bernoulli distribution
    pXgivXstarZ = expit_XgivXstarZ ^ complete_data[, X_name] * (1 - expit_XgivXstarZ) ^ (1 - complete_data[, X_name])
  }

  ## P(Y, X|X*, Z) OR P(Y,X|X*)
  pYXgivXstarZ = pYgivXZ * pXgivXstarZ

  ## Marginalize X out of P(Y, X|X*, Z) for unqueried
  marg_pYXgivXstarZ = rowsum(x = pYXgivXstarZ,
                             group = complete_data[, "id"])

  ### Dazzle fix: replace with another VERY small number that's close to 0
  pYgivXZ[which(pYgivXZ == 0)] = 5e-324
  pXgivXstarZ[which(pXgivXstarZ == 0)] = 5e-324
  marg_pYXgivXstarZ[which(marg_pYXgivXstarZ == 0)] = 5e-324

  # Compute log-likelihood
  ll = sum(log(pYgivXZ[c(1:n)])) +
    sum(log(pXgivXstarZ[c(1:n)])) +
    sum(log(marg_pYXgivXstarZ[-c(1:n)]))
  if(verbose) {print(paste("Queried:", ll))}

  ll = ll +
    sum(log(marg_pYXgivXstarZ[-c(1:n)]))
  if(verbose) {print(paste("Queried + Unqueried:", ll))}
  return(-ll) ## return (-1) x log-likelihood for maximization
}

make_complete_data = function(data, analysis_formula, error_formula, 
                              rows, Y, X,
                              offset, x = NULL) {
  if (!is.null(x)) { ## If forcing a particular value of X
    data[, X] = x
  }
  comp_dat = model.matrix(object = analysis_formula, 
                          data = data[rows, ]) ### Model matrix incl. intercept
  comp_dat = cbind(comp_dat, 
                   model.matrix(object = error_formula, 
                                data = data[rows, ])) ### Model matrix incl. intercept
  comp_dat = cbind(comp_dat, data[rows, c(Y, offset, "row_num")]) ### Bring in (Y, Offset, row_num)
  comp_dat = comp_dat[, unique(colnames(comp_dat))] ### Get rid of potential duplicate columns
  return(comp_dat)
}