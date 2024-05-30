
#####################################################################
#                                                                   #
#     This script contains the algorithm for running the            #
#     ML fisher-scoring algorithm in a block-multivariate           #
#     normal gaussian where the covariance matrix is block-         #
#     diagonal and a linear combination of known semi-              #
#     definite matrices                                             #
#                                                                   #
#                                                                   #
#     Dependicies: To run this script you have to source            #
#     function_library.R which contains loading of rele-            #
#     vant packages and relevant helper-functions.                  #
#                                                                   #
#####################################################################



#-------------------------------------------
#       ML fisher scoring function
#-------------------------------------------
score_fisher_function <- function(design_matrix, semi_def_matrix, outcomes, params, small_value_threshold = 1e-12, add_small_constant = 1e-9){
  
  #-------------Parameters-----------------
  beta <- params[1:ncol(design_matrix)]                             #extracting mean-value parameters
  sigma2_vec <- params[(ncol(design_matrix) + 1):length(params)]    #extracting variance parameters
  
  
  # Calculating omega inverse
  omega <- omega_func(semi_def_matrix, sigma2_vec)
  
  # Setting very small values to 0
  omega[omega < small_value_threshold] <- 0
  
  # Adding small value to diagonal if diagonal values are very small
  omega <- omega + (diag(omega) < small_value_threshold) * add_small_constant * diag(length(diag(omega)))
  
  # Inverting omega
  omega_inv <- chol2inv(chol(omega))
  
  
  
  #--------------Calculating Py matrix-----------------
  Py <- Py_func(omega_inv = omega_inv, outcomes = outcomes, design_matrix = design_matrix, beta = beta)
  
  
  
  #-------------Calculating mean value parameters----- (27.21)
  M <- t(design_matrix) %*% omega_inv %*% design_matrix
  
  
  
  #-------------Calculating S matrix-----------------
  S <- S_matrix_function(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv)
  
  
  #-------------Calculating scores-------------------
  
  
  score <- parameter_score(design_matrix = design_matrix, 
                           omega_inv = omega_inv,
                           semi_def_matrix = semi_def_matrix,
                           Py = Py, 
                           outcomes = outcomes, 
                           beta = beta, 
                           sigma2 = sigma2)
  
  
  return(list('M' = M, 'S' = S, 'score' = score))
}




#-------------------------------------------
#       ML fisher scoring algorithm
#-------------------------------------------
find_mle_parameters <- function(init_params, design_matrices, semi_def_matrices, outcome_list, update_step_size = 1, max_iter = 10000, tolerance = 1e-12, 
                                small_value_threshold = 1e-12, add_small_constant = 1e-12){
  
  
  for (iter in 1:max_iter) {
    
    out <- Map(score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_params))
    
    
    # Sum blocks
    M_sum <- Reduce('+',lapply(out, function(x) x$M))
    S_sum <- Reduce('+',lapply(out, function(x) x$S))
    
    
    # Setting very small values to 0
    M_sum[M_sum < small_value_threshold] <- 0
    S_sum[S_sum < small_value_threshold] <- 0
    
    # Adding small value to diagonal if diagonal values are very small
    M_sum <- M_sum + (diag(M_sum) < small_value_threshold) * add_small_constant * diag(length(diag(M_sum)))
    S_sum <- S_sum + (diag(S_sum) < small_value_threshold) * add_small_constant * diag(length(diag(S_sum)))
    
    
    # Define inverse fisher information
    fisher_inv <- bdiag(chol2inv(chol(M_sum)), chol2inv(chol(S_sum)))
    
    # Sum scores
    score <- rowSums(sapply(out, function(x) x$score))
    
    #Calculate update step
    update_step <- fisher_inv %*% score
    
    
    # Check for convergence
    if (sum((update_step)^2) < tolerance) {
      break
    }
    
    # Update parameters for the next iteration
    init_params <- init_params + update_step_size * update_step
  }
  
  return(init_params)
}

