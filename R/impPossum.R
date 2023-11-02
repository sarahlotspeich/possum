#' Multiple imputation for Poisson regression problems with covariate measurement error
#' This function returns the multiple imputation (MI) estimates for the Poisson regression model with covariate measurement error from Lotspeich et al. (2023+)
#'
#' @param imputation_formula imputation model formula (or coercible to formula) passed through to \code{lm()}, a formula expression as for other regression models. The response should be the error-prone version of the covariate. 
#' @param analysis_formula analysis model formula (or coercible to formula) passed through to \code{glm()}, a formula expression as for other regression models. The response should be the Poisson model outcome, and, if needed, the offset can be included in this formula using the \code{offset()} function.
#' @param data dataset containing at least the variables included in \code{imputation_formula} and \code{analysis_formula}.
#' @param B desired number of imputations. Default is \code{B = 1}, which is single imputation.
#' @param seed (optional) random seed to use for reproducibility of random draws. Default is \code{seed = NULL}, which does not reset the random seed inside the function. 
#' @return dataframe with final coefficient and standard error estimates for the analysis model, pooled according to Rubin's rules.
#' @export

impPossum = function(imputation_formula, analysis_formula, data, B = 1, seed = NULL) {
  # Fit complete case model (Poisson)
  ## Based on user-supplied analysis_formula
  ## To get the dimension of the analysis model parameters 
  cc_fit =  glm(formula = analysis_formula,
                family = poisson,
                data = data)
  dim_beta = length(coefficients(cc_fit))

  # Fit imputation model (linear)
  ## Based on user-supplied imputation_formula
  imp_mod = lm(formula = as.formula(imputation_formula), 
               data = data)
  
  ## Extract means for conditional distribution from imputation model
  mu = predict(object = imp_mod, 
               newdata = data) 
  
  ## Save the name of the variable being imputed
  imp_var = gsub("~.*", "", as.character(imputation_formula))[2]
  
  ## Extract means for conditional distribution from imputation model
  mu = predict(object = imp_mod, 
               newdata = data) 
  
  # Reorder rows of data to have non-missing values first 
  data = data[order(is.na(data[, imp_var])), ]
  N = nrow(data) ## total number of rows (missing or non-missing values)
  n = sum(!is.na(data[, imp_var])) ## number of rows with non-missing values, now ordered to be first 
  
  # Multiple imputation
  ## Build matrix to hold coefficient estimates and variances from each imputation
  imp_params = imp_vars = matrix(data = NA, 
                                 nrow = B, 
                                 ncol = dim_beta)
  
  ## If user specified random seed, reset
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ## Loop over the B iterations of imputation
  for (b in 1:B) {
    ### Draw imputed values from distribution (based on imputation model) for rows with missing values
    data[-c(1:n), imp_var] = rnorm(n = (N - n), 
                                   mean = mu[-c(1:n)], 
                                   sd = sigma(imp_mod))
    
    ### Fit outcome model with imputed X (Poisson)
    fit =  glm(formula = analysis_formula, 
               family = poisson,
               data = data)
    
    ### Save parameters
    imp_params[b, ] = coefficients(fit)
    
    ### Save standard errors
    imp_vars[b, ] = diag(vcov(fit))
  }
  
  # Save pooled parameter estimates
  res = data.frame(Coefficient = names(coefficients(fit)), 
                   Estimate = colMeans(imp_params), 
                   Standard.Error = sapply(X = 1:dim_beta, 
                                           FUN = function(c) sqrt(mean(imp_vars[, c]) + (B + 1) * mean((imp_params[, c] - mean(imp_params[, c])) ^ 2))))
  
  # Return results 
  return(res)
}