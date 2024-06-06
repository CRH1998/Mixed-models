


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




#Simulate simple data with mean value 0 and calculate RML
n_clusters = 1
n_individuals_in_cluster = 100

#Generate large dataset
family_dataset <- family_dataset_generator(n_clusters = n_clusters, n_individuals_in_cluster = n_individuals_in_cluster,
                                           n_mean_param = 1, n_variance_param = 2, mean_val = 0,
                                           sigma_0 = 3, sigma_1 = 4,
                                           seed = NA)


design_matrices <- family_dataset$design_matrices
semi_def_matrices <- family_dataset$semi_def_matrices
outcome <- family_dataset$outcome_list

X <- design_matrices[[1]]
semi_def_matrix <- semi_def_matrices[[1]]
y <- outcome[[1]]

#Calculating restricted log-likelihood
REML_ll <- function(X, semi_def_matrix, y, params){
  
  
  n_i <- nrow(X)
  m <- ncol(X)
  
  #Calculating inverse covariance matrix
  V <- omega_func(semi_def_matrix, params)
  
  #Inverse omega
  Vinv <- chol2inv(chol(V))
  
  
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X
  XtVinvX_inv <- solve(XtVinvX)
  
  P <- Vinv - Vinv %*% X %*% XtVinvX_inv %*% XtVinv
  
  yPy <- t(y) %*% P %*% y
  
  log_det_V <- sum(log(eigen(V, symmetric = T)$values))
  log_det_XtVinvX <- sum(log(eigen(XtVinvX_inv, symmetric = T)$values))
  
  
  restricted_logLik <- -0.5*(log_det_V + log_det_XtVinvX + yPy + (n_i-m)*log(2*pi))
  
  return(restricted_logLik)
}

REML_ll(X, semi_def_matrix, y, n_i, m, c(3.249801,0))

lower_bounds <- c(0.000001, 0.000001)
upper_bounds <- c(Inf, Inf)

optim(par = c(2,2), 
      fn = REML_ll, 
      X = X, semi_def_matrix = semi_def_matrix, y = y, m = m, n_i = n_i,
      lower = lower_bounds,
      upper = upper_bounds,
      method = 'L-BFGS-B',
      control=list(fnscale=-1))










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
  res <- Map(restricted_log_likelihood_block, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(parameters))
  
  return(Reduce('+', res))
}






