


########################################################
#                                                      #
# Functions for calculating restricted log-likelihood  #
# in a mixed-model with block-diagonal covariance      #
# matrix                                               #
#                                                      #
#                                                      #
# Dependencies:                                        #
# Source function_library.R                            #
#                                                      #
########################################################


#Calculating restricted log-likelihood
REML_ll <- function(X, semi_def_matrix, y, params){
  
  n <- nrow(X)
  p <- ncol(X)
  
  beta <- params[1:m]
  sigma2 <- params[m+1, length(params)]
  
  #Calculating inverse covariance matrix
  V <- omega_func(semi_def_matrix, sigma2)
  
  #Inverse omega
  Vinv <- chol2inv(chol(V))
  
  #Log determinant V
  log_det_V <- log(det(V))
  
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X
  log_det_XtVinvX <- log(det(XtVinvX))

  #ytVy
  RSS <- t(y - X %*% beta) %*% Vinv %*% (y - X %*% beta)
  log_ytVinvy <- log(t(y - X %*% beta) %*% Vinv %*% (y - X %*% beta))
  
  #log_det_V <- sum(log(eigen(V, symmetric = T)$values))
  #log_det_XtVinvX <- sum(log(eigen(XtVinvX, symmetric = T)$values))
  
  lnum <- log(2 * pi * RSS)
  nmp <- n-p
  ldw <- n
  ldL2 <- 
  
  restricted_logLik <- -0.5*((n_i-m)*log(2*pi) + log_det_V + log_det_XtVinvX + ytVinvy)
  
  return(restricted_logLik)
}





#Calculating restricted log-likelihood
REML_p_ll_blocks <- function(X, semi_def_matrix, y, params){
  
  n_i <- nrow(X)
  m <- ncol(X)
  
  #Calculating inverse covariance matrix
  V <- omega_func(semi_def_matrix, params^2)
  
  #Log determinant V
  log_det_V <- log(det(V))

  
  #Inverse omega
  Vinv <- chol2inv(chol(V))
    
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X
  
  #ytVy
  ytVinvy <- t(y) %*% Vinv %*% y
  
  
  return(list('ytVinvy' = ytVinvy, 'XtVinvX' = XtVinvX, 'log_det_V' = log_det_V))
}






#----------------------------------------------
#       Block wise restricted log likelihood function
#----------------------------------------------


restricted_log_likelihood_block <- function(design_matrix, semi_def_matrix, outcomes, parameters){
  
  #Determining n_i
  n_i <- nrow(design_matrix)
  
  sigma2_vec <- parameters
  
  
  #Calculating inverse covariance matrix
  omega <- omega_func(semi_def_matrix, sigma2_vec)
  
  #Inverse omega
  omega_inv <- chol2inv(chol(omega))
  
  #Storing repeated matrix multiplications
  XVinvX <- t(design_matrix) %*% omega_inv %*% design_matrix
  
  
  # Calculate beta_hat
  beta_hat <- chol2inv(chol(XVinvX)) %*% t(design_matrix) %*% omega_inv %*% outcomes
  yXB <- outcomes - design_matrix %*% beta_hat
  
  #Calculating log-likelihood
  res <- -0.5 * ((n_i - length(beta_hat)) * log(2 * pi) + log(det(omega)) + log(det(XVinvX)) + t(yXB) %*% omega_inv %*% (yXB))
    
  
  #Returning log-likelihood
  return(res)
}




#--------------------------------------------------------------------------------------------
#       Calculate full log-likelihood using block-wise log-likelihood
#--------------------------------------------------------------------------------------------

restricted_log_likelihood <- function(design_matrices, semi_def_matrices, outcome_list, parameters){
  
  #Applying log-likehood function to each element of lists, parameter vector and sigma vector is the same for each individual
  res <- Map(REML_ll, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(parameters))
  
  return(Reduce('+', res))
}






