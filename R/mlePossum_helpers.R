#score function (DOES NOT HANDLE NULL Z APPROPRIATELY)
U = function(beta, X, Y, Z = NULL, data) {
  
  #turn inputs into matrices
  beta = matrix(data = beta, ncol = 1)
  data = data.matrix(frame = data) #holds X, Z, and other potential covariates
  
  #initialize lambda, score
  if(is.null(Z)){
    l = exp(beta[1] + beta[2] * data[, X]) 
  }
  else {
    l = exp(beta[1] + data[, c(X, Z)] %*% beta[-1]) ## dimension n x 1
  }
  score = matrix(data = 0, 
                  nrow = nrow(beta)) ## dimension p x 1 (where p = dim(beta))
  
  # save Y - lambda 
  yml = data[, Y] - l ## dimension n x 1
  
  #fill first element of score
  score[1] = sum(yml)
  
  #save column copies of Y - lambda in a matrix
  ## same column, replicated once for each covariate in (X, Z)
  yml_wide = matrix(data = yml, 
                     nrow = length(yml),
                     ncol = length(c(X, Z)), 
                     byrow = FALSE) ## dimension n x (p - 1)
  
  #fill remaining elements of score (clever way)
  score[-1] = colSums(data[, c(X, Z)] * yml_wide)
  
  #return score
  score
}

#Fisher information
info = function(beta, X, Z = NULL, data) {
  #turn inputs into matrices
  beta = matrix(data = beta, ncol = 1)
  data = data.matrix(frame = data) #holds X, Z, and other potential covariates
  
  #initialize lambda
  if(is.null(Z)) {
    l = exp(beta[1] + beta[2] * data[, X])  
  } else {
    l = exp(beta[1] + data[, c(X,Z)] %*% beta[-1]) #dim n x 1
  }
  
  #save top left 2x2
  mini_top_left = sum(l)
  mini_top_right = sum(data[, X] * l)
  mini_bottom_left = mini_top_right
  mini_bottom_right = sum(data[,X]^2 * l)
  
  #create empty information matrix 
  i = matrix(data = 0, 
              nrow = nrow(beta), 
              ncol = nrow(beta))
  i[1, 1] = mini_top_left
  i[1, 2] = i[2, 1] = mini_bottom_left
  i[2, 2] = mini_bottom_right 
  
  if (!is.null(Z)) {
    #extend bottom left vector (sum of product of Z and lambda)
    l_wide = matrix(data = l, 
                     nrow = length(l),
                     ncol = length(Z), 
                     byrow = FALSE) 
    bottom_left = colSums(data[, Z] * l_wide)
    i[-c(1:2), 1] = bottom_left
    i[1, -c(1:2)] = t(bottom_left)
    
    #extend bottom middle vector (sum of product of X, Z, and lambda)
    xl = data[,X] * l
    xl_wide = matrix(data = xl,
                      nrow = length(xl),
                      ncol = length(Z),
                      byrow = FALSE)
    bottom_middle = colSums(data[, Z] * xl_wide)
    i[-c(1:2), 2] = bottom_middle
    i[2, -c(1:2)] = t(bottom_middle)
    
    #extend bottom right matrix
    #zzt = data[, Z] %*% t(data[, Z]) #issue in dimensionality
    zz = data[,Z] * data[,Z] #temp issue fix
    
    #l_wide = matrix(data = l, 
    #/nrow = length(l),
    #ncol = ncol(zz),
    #ncol = ncol(zzt), #same issue
    #byrow = FALSE) 
    #bottom_right = colSums(zzt * l_wide) #issue in dimensionality
    bottom_right = sum(zz * l) #temp issue fix
    i[-c(1:2), -c(1:2)] = bottom_right
  }
  return(-i)
}