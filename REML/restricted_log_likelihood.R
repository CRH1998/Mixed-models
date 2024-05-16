


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
  XVX <- t(design_matrix) %*% omega_inv %*% design_matrix
  
  
  # Calculate beta_hat
  beta_hat <- chol2inv(chol(XVX)) %*% t(design_matrix) %*% omega_inv %*% outcomes
  yXB <- outcomes - design_matrix %*% beta_hat
  
  #Calculating log-likelihood
  res <- -0.5 * ((n_i - length(beta_hat)) * log(2 * pi) + log(det(omega)) + log(det(XVX)) + t(yXB) %*% omega_inv %*% (yXB))
    
  
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






