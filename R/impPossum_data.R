#' Multiple imputation for Poisson regression problems with covariate measurement error
#' This function returns the multiply imputed datasets 
#'
#' @param imputation_formula imputation model formula (or coercible to formula) passed through to \code{lm()}, a formula expression as for other regression models. The response should be the error-prone version of the covariate. 
#' @param data dataset containing at least the variables included in \code{imputation_formula} and \code{analysis_formula}.
#' @param B desired number of imputations. Default is \code{B = 1}, which is single imputation.
#' @param seed (optional) random seed to use for reproducibility of random draws. Default is \code{seed = NULL}, which does not reset the random seed inside the function. 
#' @return dataframe with final coefficient and standard error estimates for the analysis model, pooled according to Rubin's rules.
#' @export

impPossum_data = function(imputation_formula, data, B = 1, seed = NULL) {
  # Fit imputation model (linear)
  ## Based on user-supplied imputation_formula
  imp_mod = lm(formula = as.formula(imputation_formula), 
               data = data)
  
  ## Save the name of the variable being imputed
  imp_var = gsub("~.*", "", as.character(imputation_formula))[2]
  
  # Reorder rows of data to have non-missing values first 
  data = data[order(is.na(data[, imp_var])), ]
  N = nrow(data) ## total number of rows (missing or non-missing values)
  n = sum(!is.na(data[, imp_var])) ## number of rows with non-missing values, now ordered to be first 
  
  ## Extract means for conditional distribution from imputation model
  mu = as.vector(predict(object = imp_mod, 
                         newdata = data)) 
  
  if (B > 0) {
    # Stochastic imputation
    ## Build list to hold imputed datasets
    imp_data = list()
    
    ## If user specified random seed, reset
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    ## Loop over the B iterations of imputation
    for (b in 1:B) {
      ### Draw imputed values from distribution for rows with missing values
      #### (based on imputation model) 
      data[-c(1:n), imp_var] = rnorm(n = (N - n), 
                                     mean = mu[-c(1:n)], 
                                     sd = sigma(imp_mod))
      
      ### Save imputed data to the list
      imp_data[[b]] = data
    }
    
    # Return results 
    return(imp_data)
  } else {
    # Deterministic imputation
    data[-c(1:n), imp_var] = mu[-c(1:n)]
    
    # Return results 
    return(data)
  }
}