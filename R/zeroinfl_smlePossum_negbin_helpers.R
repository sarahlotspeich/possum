#' @importFrom bizicount dzinb
get_pYgivX_from_fit = function(fit, beta, eta, theta, data, beta_cols, eta_cols, Y) {
  # If parameter values weren't supplied, take them from the fitted model
  data = data.matrix(data)
  if (is.null(beta)) {
    beta = matrix(data = fit$coef[grepl(pattern = "ct_", x = names(fit$coef))],
                  ncol = 1)
    eta = matrix(data = fit$coef[grepl(pattern = "zi_", x = names(fit$coef))],
                 ncol = 1)
    theta = fit$coef["Theta"]
  } else {
    beta = matrix(data = beta,
                  ncol = 1)
    eta = matrix(data = eta,
                 ncol = 1)
  }
  dzinb(
    x = data[, Y],
    size = theta,
    mu = exp(data[, beta_cols] %*% beta),
    psi =  1 / (1 + exp(-data[, eta_cols] %*% eta))
  )
}
