#' Hessian Row
#' 
#' Multiply two variables and return colsums
#' TODO
#' 
#' @param x matrix
#' @param pm post multiply
#' @noRd
hessian_row <- function(x, pm) {
  return(colSums(x * pm))
}
