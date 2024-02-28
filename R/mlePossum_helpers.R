loglik_mat = function(beta_eta, 
                      Y_name, X_name, 
                      Z_name = NULL, Xstar_name, 
                      Q_name, data,
                      verbose = FALSE) {
  #print(beta_eta)
  
  # Save useful constants
  N = nrow(data) ## Phase I sample size
  n = sum(data[, Q_name]) ## Phase II sample size
  
  # Reorder data to put queried rows first
  data = data[order(data[, Q_name], decreasing = TRUE), ]
  
  # Create matrix of complete data
  if (n < N) {
    queried_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Z_name, Xstar_name)])
    unqueried_data = rbind(
      cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 0, data[-c(1:n), c(Z_name, Xstar_name)]),
      cbind(id = (n+1):N, data[-c(1:n), Y_name], X_name = 1, data[-c(1:n), c(Z_name, Xstar_name)])
    )
    colnames(unqueried_data) = c("id", Y_name, X_name, Z_name, Xstar_name)
    complete_data = data.matrix(rbind(queried_data, unqueried_data))
  } else {
    complete_data = cbind(id = 1:n, data[1:n, c(Y_name, X_name, Z_name, Xstar_name)])
  }
  
  # Compute log-likelihood 
  ## P(Y|X,Z) from Poisson distribution
  lambdaY = exp(beta_eta[1] + beta_eta[2] * complete_data[, X_name] + beta_eta[3] * complete_data[, Z_name])
  
  ### Dazzle fix: replace y with data[, Y_name]
  pYgivXZ = dpois(x = complete_data[, Y_name], lambda = lambdaY)
  
  ## P(X|X*,Z) from Bernoulli distribution
  pXgivXstarZ = 1 / (1 + exp(-(beta_eta[4] + beta_eta[5] * complete_data[, Xstar_name] + beta_eta[6] * complete_data[, Z_name]))) ^ complete_data[, X_name] * 
    (1 - 1 / (1 + exp(-(beta_eta[4] + beta_eta[5] * complete_data[, Xstar_name] + beta_eta[6] * complete_data[, Z_name])))) ^ (1 - complete_data[, X_name]) 
  
  ## P(Y, X|X*, Z) 
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