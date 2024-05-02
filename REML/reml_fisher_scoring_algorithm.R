










#-------------------------------------------
#       Fisher scoring function
#-------------------------------------------
reml_score_fisher_function <- function(design_matrix, semi_def_matrix, outcomes, params){
  

  # Calculating omega inverse
  omega <- omega_func(semi_def_matrix = semi_def_matrix, sigma2_vec = params)
  omega_inv <- chol2inv(chol(omega))
  
  
  # Calculating P matrix
  P <- P_func(omega_inv = omega_inv, design_matrix = design_matrix)
  

  #-------------Calculating S matrix-----------------
  S <- S_matrix_reml_function(semi_def_matrix = semi_def_matrix, P = P)
  
  
  #-------------Calculating scores-------------------
  score <- reml_score_func(P = P, outcomes = outcomes, semi_def_matrix = semi_def_matrix)
  
  
  return(list('S' = S, 'score' = score))
}




#-------------------------------------------
#       Fisher scoring algorithm
#-------------------------------------------
find_remle_parameters <- function(init_params, design_matrices, semi_def_matrices, outcome_list, max_iter = 1000000, tolerance = 1e-3, update_step_size = 1){
  
  max_iter <- max_iter
  tolerance <- tolerance
  
  for (iter in 1:max_iter) {
    
    out <- Map(reml_score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_params))

    # Sum blocks
    S_sum <- Reduce('+',lapply(out, function(x) x$S))
    
    #print('S_sum')
    #print(S_sum)
    
    # Define inverse fisher information
    fisher_inv <- chol2inv(chol(S_sum))
    
    #print('fisher_inv')
    #print(fisher_inv)
    
    # Sum scores
    score <- rowSums(sapply(out, function(x) x$score))
    
    #print('score')
    #print(score)
    
    
    #Calculate update step
    update_step <- fisher_inv %*% score
    
    #print('update_step')
    #print(update_step)
    
    # Check for convergence
    if (sum((update_step)^2) < tolerance) {
      break
    }
    
    # Update parameters for the next iteration
    init_params <- init_params + update_step_size * update_step
    
    #print('init_params')
    #print(init_params)
  }
  
  return(init_params)
}




