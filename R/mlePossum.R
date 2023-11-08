#' Maximum likelihood estimation for Poisson regression problems with covariate measurement error
#' This function returns the maximum likelihood estimates (MLEs) for the Poisson regression model with covariate measurement error from Mullan et al. (2023+)
#'
#' @param start_guess numeric vector, treated as initial guess for model coefficients
#' @param x string, name of the column in \code{data} treated as the error-prone covariate 
#' @param y string, name of the column in \code{data} treated as the outcome
#' @param z string or string vector, name(s) of the column(s) in \code{data} treated as the error-free covariates. Default is \code{z = NULL} for no error-free covariates.
#' @param data dataset containing at least the variables \code{x}, \code{y}, and \code{z}.
#' @param maxtol (optional) scalar, the maximum allowed difference between guesses used to define convergence. Default is \code{maxtol = 1E-5}.
#' @param maxiter (optional) scalar, the maximum number of guesses allowed until algorithm finishes without convergence. Default is \code{maxiter = 1E3}.
#' @param verbose (optional) logical, if \code{TRUE} progress messages are displayed throughout; if \code{FALSE} (the default) no messages are displayed.
#' @return a list of two items 
#' \item{coefficients}{dataframe with final coefficient and standard error estimates for the analysis model.}
#' \item{convergence}{string message indicating convergence status}
#' @export

mlePossum = function(start_guess, x, y, z = NULL, data, maxtol = 1E-5, maxiter = 1E3, verbose = FALSE) {
  beta_curr = start_guess
  diff = maxtol + 1
  iterations = 0
  while(diff > maxtol & iterations <= maxiter){
    beta_next = beta_curr - solve(info(beta_curr, x, z, data)) %*% U(beta_curr, x, y, z, data)
    diff = max(abs(beta_next - beta_curr))
    iterations = iterations + 1
    beta_curr = beta_next
    if(verbose) {
      print(paste("beta_next:", beta_next,"\n"))
      print(paste("diff:", diff, "\n"))
      print(paste("iterations:", iterations, "\n"))
    }
  }
  if(diff > maxtol & iterations >= maxiter){
    conv_msg = "We hit the maximum number of iterations but did not converge."
  }
  else { conv_msg = "We have achieved convergence!"}
  errors = diag(solve(info(beta_curr, x, z, data)))
  return(list(coefficients = data.frame(estimates = beta_curr,
                                        std.errors = errors),
              convergence = conv_msg))
}