
#####################################################################
#                                                                   #
#     This script contains the algorithm for running the            #
#     REML fisher-scoring algorithm in a block-multivariate         #
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
#       REML fisher scoring function
#-------------------------------------------
reml_score_fisher_function <- function(design_matrix, semi_def_matrix, outcomes, params, small_value_threshold = 1e-12, add_small_constant = 1e-12){
  

  # Calculating omega inverse
  omega <- omega_func(semi_def_matrix = semi_def_matrix, sigma2_vec = params)
  
  # Setting very small values to 0
  omega[omega < small_value_threshold] <- 0
  
  # Adding small value to diagonal if diagonal values are very small
  omega <- omega + (diag(omega) < small_value_threshold) * add_small_constant * diag(length(diag(omega)))
  
  # Inverting omega
  omega_inv <- chol2inv(chol(omega))
  
  # Calculating P matrix
  P <- P_func(omega_inv = omega_inv, design_matrix = design_matrix)
  
  
  #Multiplying the P matrix with each semi-definite matrix (variance component matrix) to save computations
  A <- multiply_list_by_matrix(P, semi_def_matrix)

  
  #-------------Calculating S matrix-----------------
  S <- S_matrix_reml_function(P = P, semi_def_matrix = semi_def_matrix)
  
  
  #-------------Calculating scores-------------------
  score <- reml_score_func(P = P, outcomes = outcomes, semi_def_matrix = semi_def_matrix)
  
  
  return(list('S' = S, 'score' = score))
}




#-------------------------------------------
#       REML fisher scoring algorithm
#-------------------------------------------
find_remle_parameters <- function(init_params, design_matrices, semi_def_matrices, outcome_list, max_iter = 1000000, 
                                  tolerance = 1e-6, update_step_size = 1, small_value_threshold = 1e-12, add_small_constant = 1e-6){
  
  for (iter in 1:max_iter) {
    print(iter)
    out <- Map(reml_score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_params))

    # Sum blocks
    S_sum <- Reduce('+',lapply(out, function(x) x$S))
    
    # Setting very small values to 0
    S_sum[S_sum < small_value_threshold] <- 0
    
    # Adding small value to diagonal if diagonal values are very small
    S_sum <- S_sum + (S_sum < small_value_threshold) * add_small_constant
    
    
    # Define inverse fisher information
    fisher_inv <- chol2inv(chol(S_sum))
    
    
    # Sum scores
    score <- rowSums(sapply(out, function(x) x$score))
    
    
    #Calculate update step
    update_step <- fisher_inv %*% score
    
    print(update_step)
    
    # Check for convergence
    if (sum((update_step)^2) < tolerance) {
      break
    }
    
    # Update parameters for the next iteration
    init_params <- init_params + update_step_size * update_step
  }
  
  init_params[init_params < 0] <- 1e-12
  
  return(init_params)
}




