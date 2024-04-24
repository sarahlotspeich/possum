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
