#' score
#'
#' @param beta a numeric vector of parameters
#' @param x a numeric vector of covariates
#' @param y a numeric vector of outcomes
#'
#' @return a scalar
#' @export
#'
#score function
U <- function(beta, x, y) {
  
  #initialize lambda, score, sample size
  l <- exp(beta[1] + beta[2]*x)
  score <- rep(0, times = 2)
  
  #fill first element of score
  score[1] <- sum(y - l)
  
  #fill second element of score
  score[2] <- sum(x * (y - l))
  
  #return score
  score
}