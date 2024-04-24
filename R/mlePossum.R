#' Maximum likelihood estimation for Poisson regression problems with covariate misclassification
#' This function returns the maximum likelihood estimates (MLEs) for the Poisson regression model with covariate misclassification from Mullan et al. (2024+)
#'
#' @param error_formula misclassification model formula (or coercible to formula), a formula expression as for other regression models. The response should be the error-free version of the error-prone of the covariate.
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be the Poisson model outcome, and, if needed, the offset can be provided as the \code{offset} argument.
#' @param offset optional, variable name string for the analysis model offset. Default is \code{offset = NULL} for no offset.
#' @param data dataset containing at least the variables included in \code{error_formula} and \code{analysis_formula}.
#' @param noFN logical, if \code{noFN = FALSE} (the default), then it is assumed that there can be both false positives and false negatives in the error-prone exposure. If \code{noFN = TRUE}, the error mechanism is restricted to only false positives.
#' @return dataframe with final coefficient and standard error estimates for the analysis model
#' @export

mlePossum = function(error_formula, analysis_formula, offset = NULL, data, noFN = FALSE) {
  ## Extract variable names from user-specified formulas
  get_Y_name = as.character(as.formula(analysis_formula))[2]
  get_X_name = as.character(as.formula(error_formula))[2]

  analysis_covar = unlist(strsplit(x = gsub(pattern = " ",
                                            replacement = "",
                                            x = as.character(as.formula(analysis_formula))[3]),
                                   split = "+",
                                   fixed = TRUE))
  error_covar = unlist(strsplit(x = gsub(pattern = " ",
                                         replacement = "",
                                         x = as.character(as.formula(error_formula))[3]),
                                split = "+",
                                fixed = TRUE))
  get_Xstar_name = setdiff(error_covar, analysis_covar)

  get_Z_name = intersect(error_covar, analysis_covar)

  ## Fit complete-case models to get initial values
  ### P(Y|X,Z)
  if (is.null(offset)) {
    cc_fit = glm(formula = as.formula(analysis_formula),
                 data = data,
                 family = poisson)$coefficients
  } else {
    cc_fit = glm(formula = as.formula(paste0(paste(get_Y_name, paste(get_X_name, collapse = "+"), sep = "~"), "+offset(log(", offset, "))")),
                 data = data,
                 family = poisson)$coefficients
  }
  ### P(X|X*,Z)
  cc_fit = c(cc_fit,
             glm(formula = as.formula(error_formula),
                 data = data,
                 family = binomial)$coefficients)

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
                      offset_name = offset,
                      data = data,
                      noFN = noFN)
  } else {
    optim_res = optim(fn = loglik_mat,
                      par = cc_fit,
                      hessian = TRUE,
                      method = "BFGS",
                      Y_name = get_Y_name,
                      X_name = get_X_name,
                      Xstar_name = get_Xstar_name,
                      offset_name = offset,
                      Q_name = "Q",
                      data = data,
                      noFN = noFN)
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
  return(list(coefficients = res,
              #misclass_coefficients = misclass_res,
              convergence = optim_res$convergence))
}
