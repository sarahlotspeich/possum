#' Maximum likelihood estimation for Poisson regression problems with covariate misclassification
#' This function returns the maximum likelihood estimates (MLEs) for the Poisson regression model with covariate misclassification from Mullan et al. (2024+)
#'
#' @param error_formula misclassification model formula (or coercible to formula), a formula expression as for other regression models. The response should be the error-free version of the error-prone of the covariate. 
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be the Poisson model outcome, and, if needed, the offset can be included in this formula using the \code{offset()} function.
#' @param data dataset containing at least the variables included in \code{error_formula} and \code{analysis_formula}.
#' @return dataframe with final coefficient and standard error estimates for the analysis model
#' @export

mlePossum = function(error_formula, analysis_formula, data) {
  ## Fit complete-case models to get initial values 
  ### P(Y|X,Z) 
  cc_fit = glm(formula = as.formula(analysis_formula), 
               data = data, 
               family = poisson)$coefficients
  ### P(X|X*,Z)
  cc_fit = c(cc_fit, 
             glm(formula = as.formula(error_formula), 
                 data = data, 
                 family = binomial)$coefficients)

  ## Extract variable names from user-specified formulas 
  get_Y_name = as.character(as.formula(analysis_formula))[2]
  get_X_name = as.character(as.formula(error_formula))[2]
  
  analysis_covar = unlist(strsplit(x = gsub(pattern = " ", replacement = "", x = as.character(as.formula(analysis_formula))[3]), split = "+", fixed = TRUE))
  error_covar = unlist(strsplit(x = gsub(pattern = " ", replacement = "", x = as.character(as.formula(error_formula))[3]), split = "+", fixed = TRUE))
  get_Xstar_name = setdiff(error_covar, analysis_covar) 
  
  get_Z_name = intersect(error_covar, analysis_covar) 

  ## Add queried/non-missing data indicator
  data[, "Q"] = as.numeric(!is.na(data[, get_X_name]))
  
  if (length(get_Z_name) > 0) {
    optim_res = optim(fn = loglik_mat, 
                      par = cc_fit, 
                      hessian = TRUE, 
                      method = "BFGS",
                      Y_name = get_Y_name,
                      X_name = get_X_name,
                      Z_name = get_Z_name,
                      Xstar_name = get_Xstar_name,
                      Q_name = "Q",
                      data = data)
  } else {
    optim_res = optim(fn = loglik_mat, 
                      par = cc_fit, 
                      hessian = TRUE, 
                      method = "BFGS",
                      Y_name = get_Y_name,
                      X_name = get_X_name,
                      Xstar_name = get_Xstar_name,
                      Q_name = "Q",
                      data = data)
  }
  
  ## Prepare model output to be returned 
  ### Invert the hessian to get estimated standard errors
  num_analysis_covar = length(analysis_covar) + 1
  cov_theta = tryCatch(expr = solve(optim_res$hessian)[1:num_analysis_covar, 1:num_analysis_covar],
                       error = function(err) {
                         matrix(data = NA, 
                                nrow = num_analysis_covar, 
                                ncol = num_analysis_covar)
                       })
  
  if (optim_res$convergence == 0) {
    res = data.frame(Est = optim_res$par[1:num_analysis_covar],
                     SE = sqrt(diag(cov_theta)))
    rownames(res) = c("(Intercept)", analysis_covar)
  } else {
    res = data.frame(Est = rep(NA, num_analysis_covar),
                     SE = rep(NA, num_analysis_covar))
    rownames(res) = c("(Intercept)", analysis_covar)
  }
  
  ## Return model output
  return(res)
}