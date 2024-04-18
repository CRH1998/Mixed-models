# Loading relevant packages
library(Matrix)               # For matrix manipulation
library(clusterGeneration)    # For positive semi-definite matrices
library(mvtnorm)              # For log-likelihood calculation
library(lme4)                 # For mixed model
library(microbenchmark)       # For testing function speed
library(psych)                # For calculating the trace
library(profvis)              # For evaluating performance of code





#Helper functions
# Matrix product list by matrix
matrix_mult_specify <- function(matrix1, mult_by_right = F){
  
  if(mult_by_right == F){
    return(function(matrix2){matrix1 %*% matrix2})
  } else {
    return(function(matrix2){matrix2 %*% matrix1})
  }
}

matrix_mult_list_by_matrix <- function(matrix, list, mult_by_right = F){
  matrix_mult_matrix <- matrix_mult_specify(matrix, mult_by_right)
  return(lapply(list, matrix_mult_matrix))
}




# --------------P matrix------------------ (27.17b)
P_func <- function(omega_inv, design_matrix){
  
  A <- t(design_matrix) %*% omega_inv
  
  return(omega_inv - omega_inv %*% design_matrix %*% chol2inv(chol(A %*% design_matrix)) %*% A)
}


# If P has already been calculated
Py_func_P <- function(P, outcome){
  return(P %*% outcome)
}

y_t_P_func <- function(P, outcome){
  return(crossprod(outcome, P))
}




# Calculate S matrix
S_matrix_reml_function <- function(semi_def_matrix, P){
  
  S <- matrix(data = NA, nrow = length(semi_def_matrix), ncol = length(semi_def_matrix))
  
  for (i in 1:length(semi_def_matrix)){
    for (j in 1:length(semi_def_matrix)){
      S[i,j] <- 0.5 * tr(P %*% semi_def_matrix[[i]] %*% P %*% semi_def_matrix[[j]])
    }
  }
  return(S)
}




reml_score_func <- function(P, outcomes, semi_def_matrix){
  
  # Calculating tr(P %*% V_i)
  
  #P_semi_def <- multiply_list_by_matrix(-0.5 * P, semi_def_matrix)
  #
  #trP_semi_def <- lapply(P_semi_def, FUN = tr)
  
  
  # Calculating y^T %*% P %*% V_i %*% P %*% y

  #y_t_P_semi_def <- matrix_mult_list_by_matrix(t(outcomes) %*% P, semi_def_matrix)
  #
  #y_t_P_semi_def_Py <- matrix_mult_list_by_matrix(0.5 * P %*% outcomes, y_t_P_semi_def, mult_by_right = T)
  

    
  # Calculating score of variance components
  
  sigma_scores <- list()
  
  for (i in 1:length(semi_def_matrix)){
    sigma_scores[[i]] <- 0.5 * (-tr(P %*% semi_def_matrix[[i]]) + (t(outcomes) %*% P) %*% semi_def_matrix[[i]] %*% (P %*% outcomes))
  }
  
  #sigma_scores <- Map('+', y_t_P_semi_def_Py, trP_semi_def)
  
  return(c(unlist(sigma_scores)))
}



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
    
    print('S_sum')
    print(S_sum)
    
    # Define inverse fisher information
    fisher_inv <- chol2inv(chol(S_sum))
    
    print('fisher_inv')
    print(fisher_inv)
    
    # Sum scores
    score <- rowSums(sapply(out, function(x) x$score))
    
    print('score')
    print(score)
    
    
    #Calculate update step
    update_step <- fisher_inv %*% score
    
    print('update_step')
    print(update_step)
    
    # Check for convergence
    if (sum((update_step)^2) < tolerance) {
      break
    }
    
    # Update parameters for the next iteration
    init_params <- init_params + update_step_size * update_step
    
    print('init_params')
    print(init_params)
  }
  
  return(init_params)
}




