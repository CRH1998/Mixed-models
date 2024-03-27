# Loading relevant packages
library(Matrix)               # For matrix manipulation
library(clusterGeneration)    # For positive semi-definite matrices
library(mvtnorm)              # For log-likelihood calculation
library(lme4)                 # For mixed model
library(microbenchmark)       # For testing function speed
library(psych)                # For calculating the trace
library(profvis)              # For evaluating performance of code




#-------------------------------------------
#   Generate data for testing functions
#-------------------------------------------

# Generate function input for testing
gen_func_input <- function(n_blokke, REML = FALSE){
  
  if (n_blokke == 1){
    DF <- data.frame(y = 1:5,
                     klasse=c(1,1,1,1,1))
    
    model <- lme4::lmer(y ~ 1 + (1|klasse), data=DF, REML=REML)
    summary_model <- summary(model)
    LL_model <- logLik(model)
    parameters <- summary_model$coefficients
    
    
    #----------Gruppe 1-----------------
    model_matrix1 <- matrix(matrix(rep(1,nrow(DF)))[1:5,])
    gamma1_matrix1 <- as.matrix(diag(1, nrow = 5))
    gamma2_matrix1 <- as.matrix(bdiag(matrix(1,5,5)))
    outcome1 <- DF$y[1:5]
    
    design_matrices <- list(model_matrix1)
    semi_def_matrices <- list(list(gamma1_matrix1, gamma2_matrix1))
    
    outcome_list <- list(outcome1)
  } else if (n_blokke == 2) {
    
    DF <- data.frame(y = 1:11,
                     klasse=c(1,1,1,1,1,
                              2,2,2,2,2,2))
    
    model <- lme4::lmer(y ~ 1 + (1|klasse), data=DF, REML=REML)
    summary_model <- summary(model)
    LL_model <- logLik(model)
    parameters <- summary_model$coefficients
    
    
    #----------Gruppe 1-----------------
    model_matrix1 <- matrix(matrix(rep(1,nrow(DF)))[1:5,])
    gamma1_matrix1 <- as.matrix(diag(1, nrow = 5))
    gamma2_matrix1 <- as.matrix(bdiag(matrix(1,5,5)))
    outcome1 <- DF$y[1:5]
    
    
    #----------Gruppe 2-----------------
    model_matrix2 <- matrix(matrix(rep(1,nrow(DF)))[6:10,])
    gamma1_matrix2 <- as.matrix(diag(1, nrow = 6))
    gamma2_matrix2 <- as.matrix(bdiag(matrix(1,6,6)))
    outcome2 <- DF$y[6:11]
    
    design_matrices <- list(model_matrix1, model_matrix2)
    semi_def_matrices <- list(list(gamma1_matrix1, gamma2_matrix1),
                              list(gamma1_matrix2, gamma2_matrix2))
    
    outcome_list <- list(outcome1, outcome2)
  } else {
    DF <- data.frame(y = 1:15,
                     klasse=c(1,1,1,1,1,
                              2,2,2,2,2,
                              3,3,3,3,3))
    
    model <- lme4::lmer(y ~ 1 + (1|klasse), data=DF, REML=REML)
    summary_model <- summary(model)
    LL_model <- logLik(model)
    parameters <- summary_model$coefficients
    
    
    #----------Gruppe 1-----------------
    model_matrix1 <- matrix(matrix(rep(1,nrow(DF)))[1:5,])
    gamma1_matrix1 <- as.matrix(diag(1, nrow = 5))
    gamma2_matrix1 <- as.matrix(bdiag(matrix(1,5,5)))
    outcome1 <- DF$y[1:5]
    
    
    #----------Gruppe 2-----------------
    model_matrix2 <- matrix(matrix(rep(1,nrow(DF)))[6:10,])
    gamma1_matrix2 <- as.matrix(diag(1, nrow = 5))
    gamma2_matrix2 <- as.matrix(bdiag(matrix(1,5,5)))
    outcome2 <- DF$y[6:10]
    
    
    #----------Gruppe 3-----------------
    model_matrix3 <- matrix(matrix(rep(1,nrow(DF)))[11:15,])
    gamma1_matrix3 <- as.matrix(diag(1, nrow = 5))
    gamma2_matrix3 <- as.matrix(bdiag(matrix(1,5,5)))
    outcome3 <- DF$y[11:15]
    
    
    design_matrices <- list(model_matrix1, model_matrix2, model_matrix3)
    semi_def_matrices <- list(list(gamma1_matrix1, gamma2_matrix1),
                              list(gamma1_matrix2, gamma2_matrix2),
                              list(gamma1_matrix3, gamma2_matrix3))
    
    outcome_list <- list(outcome1, outcome2, outcome3)
  }
  
  return(list('design_matrices' = design_matrices, 'semi_def_matrices' = semi_def_matrices, 'outcome_list' = outcome_list, 'summary_model' = summary_model))
}


# Generate large data set for testing
large_dataset_generator <- function(n_clusters, n_individuals_in_cluster, seed = 1){
  
  set.seed(1)
  
  klasser <- rep(seq(1,n_clusters),each = n_individuals_in_cluster)
  subklasser <- rep(seq(1,n_clusters*2), each = n_individuals_in_cluster/2)
  y <- klasser + rnorm(n_individuals_in_cluster* n_clusters, mean = 1, sd = 1) + subklasser
  
  DF <- data.frame(y = y, klasse = klasser, subklasse = subklasser)
  
  outcome_list <- split(DF$y, rep(1:n_clusters, each = n_individuals_in_cluster, length.out = length(DF$y)))
  
  #model <- lme4::lmer(y ~ 1 + (1|klasse) + (1|subklasse), data=DF, REML=F)
  #summary_model <- summary(model)
  
  gamma0_matrix <- as.matrix(diag(1, nrow = n_individuals_in_cluster))
  gamma1_matrix <- as.matrix(bdiag(matrix(1,n_individuals_in_cluster,n_individuals_in_cluster)))
  gamma2_matrix <- as.matrix(bdiag(matrix(1, nrow = floor(n_individuals_in_cluster/2), ncol = floor(n_individuals_in_cluster/2)), 
                                   matrix(1, nrow = ceiling(n_individuals_in_cluster/2), ncol = ceiling(n_individuals_in_cluster/2))))
  gamma_list <- list(gamma0_matrix, gamma1_matrix, gamma2_matrix)
  
  semi_def_matrices <- replicate(n_clusters, gamma_list, simplify = F)
  design_matrices <- replicate(n_clusters, as.matrix(rep(1, n_individuals_in_cluster)), simplify = F)
  
  return(list('design_matrices' = design_matrices, 'semi_def_matrices' = semi_def_matrices, 'outcome_list' = outcome_list, 'DF' = DF))#, 'summary_model' = summary_model))
}


#-------------------------------------------
#   Helper functions for calculations
#-------------------------------------------

# Calculate omega
omega_func <- function(semi_def_matrix, sigma2_vec){
  
  #1) Multiplying each semi-definite matrix with sigma2
  omega_temp <- lapply(seq_along(semi_def_matrix), function(i) {
    semi_def_matrix[[i]] * sigma2_vec[i]
  })
  
  #2) Summing each term in omega_temp
  omega <- Reduce('+', omega_temp)
  
  return(omega)
}



# Multiply list by matrix
matrix_mult <- function(matrix1){
  return(function(matrix2){matrix1 %*% matrix2})
}

multiply_list_by_matrix <- function(matrix, list){
  
  matrix_mult_matrix <- matrix_mult(matrix)
  
  return(lapply(list, matrix_mult_matrix))
}



# -------------S matrix----------------- (27.22)

S_matrix_function <- function(semi_def_matrix, omega_inv){
  
  #S <- function(matrix1, matrix2, S_elements) {
  #  outer(matrix1, matrix2, function(x,y) vapply(seq_along(x), function(i) S_elements(x[[i]], y[[i]]), numeric(1)))
  #}
  #
  #return(S(A, A, S_elements))
  
  #S <- 0.5 * outer(seq_along(semi_def_matrix), seq_along(semi_def_matrix),
  #            Vectorize(function(i, j) tr(A[[i]] %*% A[[j]])))
  
  
  A <- multiply_list_by_matrix(omega_inv, semi_def_matrix)
  
  S <- matrix(data = NA, nrow = length(semi_def_matrix), ncol = length(semi_def_matrix))
  
  for (i in 1:length(semi_def_matrix)){
    for (j in 1:length(semi_def_matrix)){
      S[i,j] <- 0.5 * sum(A[[i]] * A[[j]])
    }
  }
  
  return(S)
  
}



# --------------Py matrix----------------- (27.17c)

Py_func <- function(omega_inv, outcomes, design_matrix, beta){
  return(omega_inv %*% (outcomes - design_matrix %*% beta))
}



#-------------Calculating scores------------------- (27.10)
parameter_score <- function(design_matrix, semi_def_matrix, omega_inv, Py, outcomes, beta, sigma2){
  
  score_beta <- crossprod(design_matrix, omega_inv) %*% (outcomes - design_matrix %*% beta)
  
  score_sigma <- sapply(seq_along(semi_def_matrix), function(i) {
    -0.5 * sum(omega_inv * semi_def_matrix[[i]]) + 0.5 * t(Py) %*% semi_def_matrix[[i]] %*% Py
  })
  
  score <- c(score_beta, score_sigma)
  
  return(score)
}




#-------------------------------------------
#       Log likelihood functions
#-------------------------------------------


log_likelihood_block <- function(design_matrix, semi_def_matrix, outcomes, parameters){
  
  #Determining n_i
  n_i <- nrow(design_matrix)
  
  mean_vec <- parameters[1:ncol(design_matrix)]
  sigma2_vec <- parameters[(ncol(design_matrix) + 1):length(parameters)]
  
  
  #Calculating inverse covariance matrix
  omega <- omega_func(semi_def_matrix, sigma2_vec)
  
  #Inverse omega
  omega_inv <- chol2inv(chol(omega))
  
  
  #Calculating log-likelihood
  res <- -n_i/2 * log(2 * pi) - 1/2 * log(det(omega)) - 1/2 * t(outcomes - design_matrix %*% mean_vec) %*% omega_inv %*% (outcomes - design_matrix %*% mean_vec)
  
  #Returning log-likelihood
  return(res)
}



log_likelihood <- function(design_matrices, semi_def_matrices, outcome_list, parameters){
  
  #Applying log-likehood function to each element of lists, parameter vector and sigma vector is the same for each individual
  res <- Map(log_likelihood_block, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(parameters))
  
  return(Reduce('+', res))
}







#-------------------------------------------
#       Fisher scoring function
#-------------------------------------------
score_fisher_function <- function(design_matrix, semi_def_matrix, outcomes, params){
  
  #-------------Parameters-----------------
  beta <- params[1:ncol(design_matrix)]
  sigma2_vec <- params[(ncol(design_matrix) + 1):length(params)]
  

  # Calculating omega inverse
  omega <- omega_func(semi_def_matrix, sigma2_vec)
  omega_inv <- chol2inv(chol(omega))
  
  
  
  #--------------Calculating Py matrix-----------------
  Py <- Py_func(omega_inv = omega_inv, outcomes = outcomes, design_matrix = design_matrix, beta = beta)
  
  
  
  #-------------Calculating mean value parameters----- (27.21)
  M <- crossprod(design_matrix, omega_inv) %*% design_matrix
  
  
  
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
#       Fisher scoring algorithm
#-------------------------------------------
find_mle_parameters <- function(init_params, design_matrices, semi_def_matrices, outcome_list, max_iter = 1000000, tolerance = 1e-12){
  
  max_iter <- max_iter
  tolerance <- tolerance
  
  for (iter in 1:max_iter) {
    print(iter)
    
    
    out <- Map(score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_params))
    
    # Sum blocks
    M_sum <- Reduce('+',lapply(out, function(x) x$M))
    S_sum <- Reduce('+',lapply(out, function(x) x$S))
    
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
    init_params <- init_params + update_step
  }
  
  return(init_params)
}












##Matrix multiply matrix list with itself
#matrix_mult_with_self <- function(matrix_list){
#  multiply_matrices <- function(mat) {
#    mat %*% t(mat)
#  }
#  
#  # Applying function to each matrix in the list
#  self_mult <- lapply(matrix_list, multiply_matrices)
#  
#  # Print the result
#  return(self_mult)
#}







