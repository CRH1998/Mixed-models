########################## Loading relevant packages#########################

library(Matrix)               # For matrix manipulation
library(clusterGeneration)    # For positive semi-definite matrices

library(mvtnorm)              # For log-likelihood calculation
library(lme4)                 # For mixed model

library(psych)                # For calculating the trace
library(microbenchmark)       # For testing function speed
library(profvis)              # For evaluating performance of code

library(MASS)                 # For mvrnorm()
library(kinship2)             # For family data
library(corpcor)              # For pseudoinverse

library(dplyr)                # For maintenance

############################################################################# 


#####################################################################
#                                                                   #
#     Helper functions to calculate log likelihood                  #
#     for block-multivariate normal gaussian where                  #
#     the covariance matrix is a linear combination                 #
#     of known semi-definite matrices.                              #
#                                                                   #
#     This script also contains helper functions for                #
#     running ML-fisher scoring algorihtm and REML-                 #
#     fisher scoring algorithm                                      #
#                                                                   #
#                                                                   #
#     Dependicies: This script does not depend on other             #
#     scripts. On the contrary you will need to source              #
#     this script to run the other functions in this                #
#     project                                                       #
#                                                                   #
#####################################################################



#-------------------------------------------
#   Matrix multiplication helper functions
#-------------------------------------------


# Multiply list by matrix
multiply_list_by_matrix <- function(matrix, list){
  
  matrix_mult <- function(matrix1){
    return(function(matrix2){matrix1 %*% matrix2})
  }
  
  matrix_mult_matrix <- matrix_mult(matrix)
  
  return(lapply(list, matrix_mult_matrix))
}



# Multiply list by matrix and take trace of each matrix product
tr_multiply_list_by_matrix <- function(matrix, list){
  
  tr_matrix_mult <- function(matrix1){
    return(function(matrix2){sum(matrix1 * matrix2)})
  }
  
  tr_matrix_mult_matrix <- tr_matrix_mult(matrix)
  
  return(lapply(list, tr_matrix_mult_matrix))
}



# Multiply list by matrix and specify if multiplication is on the right or the left side
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





#-------------------------------------------
#         List helper functions
#-------------------------------------------

get_outcome_mean <- function(outcome){
  return(mean(as.vector(sapply(outcome, function(x){return(x)}))))
}

get_outcome_variance <- function(outcome){
  return(var(as.vector(sapply(outcome, function(x){return(x)}))))
}






#-------------------------------------------
#   ML and REML specific helper functions
#-------------------------------------------

# Calculate omega
omega_func <- function(semi_def_matrix, sigma2_vec){
  
  omega <- 0
  
  for (i in 1:length(semi_def_matrix)){
    omega <- omega + semi_def_matrix[[i]] * sigma2_vec[i]
  }
  
  return(omega)
}



# Calculate RSS
RSS_func <-  function(X, semi_def_matrix, y, params){
  
  #Calculating inverse covariance matrix
  V <- omega_func(semi_def_matrix, params)
  Vinv <- chol2inv(chol(V))
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X
  
  #Calculating betahat
  betahat <- chol2inv(chol(XtVinvX)) %*% (XtVinv %*% y)
  
  res <- y - X %*% betahat
  
  #RSS
  RSS <- t(res) %*% Vinv %*% res
  
  return(RSS)
}



# Calculate XtVinvX
XtVinvX_func <- function(X, semi_def_matrix, y, params){
  
  #Calculating inverse covariance matrix
  V <- omega_func(semi_def_matrix, params)
  Vinv <- chol2inv(chol(V))
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X
  
  return(XtVinvX)
}





##########################################################
#             ML specific helper functions               #
##########################################################


# ---------------Calculate residual-------------------
residual_function <- function(outcomes, design_matrix, beta){
  return(outcomes - design_matrix %*% beta)
}


# -------------Calculate ML S matrix----------------- (27.22)

S_matrix_function <- function(semi_def_matrix, omega_inv){
  
  # Multiply each covariance component matrix by the inverse of the covariance matrix
  A <- multiply_list_by_matrix(omega_inv, semi_def_matrix)
  
  # Construct empty matrix to fill with values 
  S <- matrix(data = NA, nrow = length(semi_def_matrix), ncol = length(semi_def_matrix))
  
  for (i in 1:length(semi_def_matrix)){
    for (j in i:length(semi_def_matrix)){
      S[i,j] <- 0.5 * sum(A[[i]] * A[[j]])
      S[j,i] <- S[i,j]
    }
  }
  
  return(S)
  
}


# --------------Py matrix----------------- (27.17c)

Py_func <- function(omega_inv, residual_vec){
  return(omega_inv %*% residual_vec)
}


#-------------Calculate ML scores------------------- 
parameter_score <- function(XtVinv, semi_def_matrix, omega_inv, Py, residual_vec){
  
  #Fixed effects score (27.10)
  score_beta <- XtVinv %*% residual_vec
  
  
  #Random effects score (27.14b, 27.17c)
  score_sigma <- rep(NA, length(semi_def_matrix))
  
  for (i in 1:length(semi_def_matrix)){
    trVinvVi <- sum(omega_inv * semi_def_matrix[[i]])
    PytViPy <- crossprod(Py, semi_def_matrix[[i]]) %*% Py
    
    score_sigma[i] <- 0.5 * (PytViPy-trVinvVi)
  }
  
  score <- c(score_beta, score_sigma)
  
  return(score)
}






##########################################################
#             REML specific helper functions             #
##########################################################

# --------------Calculate P matrix------------------ (27.17b)
P_func <- function(omega_inv, design_matrix){
  
  A <- t(design_matrix) %*% omega_inv
  
  P <- omega_inv - t(A) %*% chol2inv(chol(A %*% design_matrix)) %*% A
  
  #P <- omega_inv - omega_inv %*% design_matrix %*% solve((t(design_matrix) %*% omega_inv %*% design_matrix)) %*% t(design_matrix) %*% omega_inv
  
  return(P)
}


# If the P matrix has already been calculated, we calculate Py and y^tP using the calculated P
Py_func_P <- function(P, outcome){
  
  Py <- P %*% outcome
  
  return(Py)
}

y_t_P_func <- function(P, outcome){
  
  ytP <- crossprod(outcome, P)
  
  return(ytP)
}




# ------------- Calculate REML S matrix----------------- (27.31a) (Note that for REML the S matrix corresponds to the fisher information matrix)
S_matrix_reml_function <- function(P, semi_def_matrix){
  
  #Here the A matrix is A <- multiply_list_by_matrix(P, semi_def_matrix), multiplying the P matrix with each semi-definite matrix (variance component matrix)
  
  # Creating empty matrix to store S matrix (fisher information matrix) entries
  S <- matrix(data = NA, nrow = length(semi_def_matrix), ncol = length(semi_def_matrix))
  
  for (i in 1:length(semi_def_matrix)){
    for (j in 1:length(semi_def_matrix)){
      S[j,i] <- 0.5 * tr(P %*% semi_def_matrix[[i]] %*% P %*% semi_def_matrix[[j]])
    }
  }
  
  return(S)
}



#-------------Calculating REML scores------------------- (27.33)
reml_score_func <- function(P, outcomes, semi_def_matrix){
  
  #Here the A matrix is A <- multiply_list_by_matrix(P, semi_def_matrix), multiplying the P matrix with each semi-definite matrix (variance component matrix)
  
  # Calculating score of variance components
  sigma_scores <- c()
  
  for (i in 1:length(semi_def_matrix)){
    sigma_scores[i] <- - 0.5 * (tr(P %*% semi_def_matrix[[i]]) - (t(outcomes) %*% P %*% semi_def_matrix[[i]]) %*% (P %*% outcomes))
  }

  return(sigma_scores)
}







