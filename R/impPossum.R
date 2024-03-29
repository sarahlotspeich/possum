#' Multiple imputation for Poisson regression problems with covariate measurement error
#' This function returns the multiple imputation (MI) estimates for the Poisson regression model with covariate measurement error from Lotspeich et al. (2023+)
#'
#' @param imputation_formula imputation model formula (or coercible to formula) passed through to \code{lm()}, a formula expression as for other regression models. The response should be the error-prone version of the covariate. 
#' @param analysis_formula analysis model formula (or coercible to formula) passed through to \code{glm()}, a formula expression as for other regression models. The response should be the Poisson model outcome, and, if needed, the offset can be included in this formula using the \code{offset()} function.
#' @param data dataset containing at least the variables included in \code{imputation_formula} and \code{analysis_formula}.
#' @param adjMatrix (optional) adjacency matrix for observations in \code{data} to be passed through to \code{spaMM::fitme()}. If \code{adjMatrix} is not \code{NULL} (the default), then a spatial random effect should be included in \code{analysis_formula}.
#' @param B desired number of imputations. Default is \code{B = 1}, which is single imputation.
#' @param truncNorm (optional) distribution to be used to draw imputed values. If \code{FALSE} (the default), a normal distribution is used; if \code{TRUE}, a truncated normal distribution is used.
#' @param seed (optional) random seed to use for reproducibility of random draws. Default is \code{seed = NULL}, which does not reset the random seed inside the function. 
#' @return dataframe with final coefficient and standard error estimates for the analysis model, pooled according to Rubin's rules.
#' @export
#' @importFrom spaMM fitme
#' @importFrom truncnorm rtruncnorm

impPossum = function(imputation_formula, analysis_formula, data, adjMatrix = NULL, B = 1, truncNorm = FALSE, seed = NULL) {
  # Fit complete case model (Poisson)
  ## Based on user-supplied analysis_formula
  ## To get the dimension of the analysis model parameters 
  if (is.null(adjMatrix)) {
    cc_fit =  glm(formula = analysis_formula,
                  family = poisson,
                  data = data)
    dim_beta = length(coefficients(cc_fit))
  } else {
    cc_fit = fitme(formula = analysis_formula, 
                   family = poisson, 
                   adjMatrix = adj_matrix,
                   data = data)
    dim_beta = nrow(as.data.frame(summary(cc_fit)["beta_table"]))
  }

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
    ### Draw imputed values from distribution for rows with missing values
    #### (based on imputation model and truncNorm argument) 
    if (truncNorm) {
      data[-c(1:n), imp_var] = rtruncnorm(n = (N - n), 
                                          a = 0, 
                                          b = Inf, 
                                          mean = mu[-c(1:n)], 
                                          sd = sigma(imp_mod))
    } else {
      data[-c(1:n), imp_var] = rnorm(n = (N - n), 
                                     mean = mu[-c(1:n)], 
                                     sd = sigma(imp_mod))
    }
    
    if (is.null(adjMatrix)) {
      ### Fit outcome model with imputed X (Poisson)
      fit = glm(formula = analysis_formula, 
                family = poisson,
                data = data)
      
      ### Save parameters
      imp_params[b, ] = coefficients(fit)
      
      ### Save standard errors
      imp_vars[b, ] = diag(vcov(fit))
    } else {
      ### Fit outcome model with imputed X (Poisson)
      fit = fitme(formula = analysis_formula, 
                  family = poisson, 
                  adjMatrix = adj_matrix,
                  data = data)
 
      ### Save parameters
      imp_params[b, ] = as.data.frame(summary(fit)["beta_table"])[, 1]
      
      ### Save standard errors
      imp_vars[b, ] = as.data.frame(summary(fit)["beta_table"])[, 2] ^ 2
    }
  }
  
  # Save pooled parameter estimates
  if (is.null(adjMatrix)) { 
    res = data.frame(Coefficient = names(coefficients(fit)), 
                     Estimate = colMeans(imp_params), 
                     Standard.Error = sapply(X = 1:dim_beta, 
                                             FUN = function(c) sqrt(mean(imp_vars[, c]) + (B + 1) * mean((imp_params[, c] - mean(imp_params[, c])) ^ 2))))
  } else {
    res = data.frame(Coefficient = rownames(as.data.frame(summary(fit)["beta_table"])), 
                     Estimate = colMeans(imp_params), 
                     Standard.Error = sapply(X = 1:dim_beta, 
                                             FUN = function(c) sqrt(mean(imp_vars[, c]) + (B + 1) * mean((imp_params[, c] - mean(imp_params[, c])) ^ 2))))
  }
  
  # Return results 
  return(res)
}