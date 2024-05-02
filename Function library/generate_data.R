######################################################
#                                                    #
# Construct R-function to calculate log likelihood   #
# for block-multivariate normal gaussian where the   #
# covariance matrix is a linear combination of kno-  #
# wn semi-definite matrices                          #
#                                                    #
#                                                    #
######################################################

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
large_dataset_generator <- function(n_clusters, n_individuals_in_cluster, sigma_0 = 1, sigma_1 = 1, seed = 1){
  
  gamma0_matrix <- as.matrix(diag(1, nrow = n_individuals_in_cluster))
  gamma1_matrix <- as.matrix(bdiag(matrix(1, nrow = floor(n_individuals_in_cluster/2), ncol = floor(n_individuals_in_cluster/2)), 
                                   matrix(1, nrow = ceiling(n_individuals_in_cluster/2), ncol = ceiling(n_individuals_in_cluster/2))))
  gamma_list <- list(gamma0_matrix, gamma1_matrix)
  
  semi_def_matrices <- replicate(n_clusters, gamma_list, simplify = F)
  design_matrices <- replicate(n_clusters, as.matrix(rep(1, n_individuals_in_cluster)), simplify = F)
  
  
  set.seed(1)
  
  klasser <- rep(seq(1,n_clusters),each = n_individuals_in_cluster)
  subklasser <- rep(seq(1:(n_clusters*2)), each = n_individuals_in_cluster/2)
  
  #rep(c(rep(1, floor(n_individuals_in_cluster/2)), rep(2, ceiling(n_individuals_in_cluster/2))), each = n_clusters*2)
  
  mu_vec <- rep(5, n_individuals_in_cluster)
  
  y <- c(replicate(n_clusters, 
                   mvrnorm(n = 1, mu = mu_vec, Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1)))
  
  DF <- data.frame(y = y, klasse = klasser, subklasse = subklasser)
  
  outcome_list <- split(DF$y, rep(1:n_clusters, each = n_individuals_in_cluster, length.out = length(DF$y)))
  
  
  
  return(list('design_matrices' = design_matrices, 'semi_def_matrices' = semi_def_matrices, 'outcome_list' = outcome_list, 'DF' = DF))#, 'summary_model' = summary_model))
}