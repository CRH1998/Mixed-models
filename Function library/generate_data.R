######################################################
#                                                    #
# This code construct a test dataset for testing     #
# the fisher scoring algorithms                      #
#                                                    #
######################################################

#-------------------------------------------
#   Generate data for testing functions
#-------------------------------------------


# Generate data set for testing
dataset_generator <- function(n_clusters, n_individuals_in_cluster, n_mean_param = 3, n_cov_param = 3,
                              mean_val = 1, beta_1 = 3, beta_2 = 3, sigma_0 = 1, sigma_1 = 1, sigma_2 = 1, sigma_3 = 1, seed = 1){
  
  
  
  # Set seed if specified for sampling later on
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  
  
  #-------------------------------------------Generate covariance structure-------------------------------------------
  
  #Diagonal matrix with ones in the diagonal and dimension corresponding to the number of individuals in each cluster (intercept)
  gamma0_matrix <- as.matrix(diag(1, nrow = n_individuals_in_cluster))
  
  
  #Block diagonal matrix with two blocks of ones and dimension corresponding to the number of individuals in each cluster
  gamma1_matrix <- as.matrix(bdiag(matrix(1, nrow = floor(n_individuals_in_cluster/2), ncol = floor(n_individuals_in_cluster/2)), 
                                   matrix(1, nrow = ceiling(n_individuals_in_cluster/2), ncol = ceiling(n_individuals_in_cluster/2))))
  
    #Corresponding feature to store in dataframe
    klasser <- rep(seq(1:(n_clusters*2)), each = n_individuals_in_cluster/2)
  
  
  
  #Block diagonal matrix with four blocks of ones and dimensions corresponding to the number of individuals in each cluster
  gamma2_matrix <- as.matrix(bdiag(matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4),
                                   matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4),
                                   matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4),
                                   matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4)))
  
    #Corresponding feature to store in dataframe
    subklasser <- rep(seq(1:(n_clusters*4)), each = n_individuals_in_cluster/4)
  
    
    
  #Block diagonal matrix with four blocks of ones and dimensions corresponding to the number of individuals in each cluster
  gamma3_matrix <- as.matrix(bdiag(matrix(1, nrow = n_individuals_in_cluster/3, ncol = n_individuals_in_cluster/3),
                                   matrix(1, nrow = n_individuals_in_cluster/3, ncol = n_individuals_in_cluster/3),
                                   matrix(1, nrow = n_individuals_in_cluster/3, ncol = n_individuals_in_cluster/3)))
    
    #Corresponding feature to store in dataframe
    other <- rep(seq(1:(n_clusters*3)), each = n_individuals_in_cluster/3)
  
    
  # Collecting gamma matrices into list
  if (n_cov_param == 1){
    gamma_list <- list(gamma0_matrix)
  } else if (n_cov_param == 2){
    gamma_list <- list(gamma0_matrix, gamma1_matrix)
  } else if (n_cov_param == 3) {
    gamma_list <- list(gamma0_matrix, gamma1_matrix, gamma2_matrix)
  } else {
    gamma_list <- list(gamma0_matrix, gamma1_matrix, gamma2_matrix, gamma3_matrix)
  }
  
  
  # Defining semi_def_matrices which is a list containing n_clusters replicas of the gamma_list
  semi_def_matrices <- replicate(n_clusters, gamma_list, simplify = F)
  
  #---------------------------------------------------------------------------------------------------------------------------------
  
  
  
  
  #-------------------------------------------Generate mean value structure-------------------------------------------
  # Defining the design_matrices 
  if (n_mean_param == 1){
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals_in_cluster))), simplify = F)
    #In this case the design matrix is simply an intercept
    x_1 <- NA
    x_2 <- NA
    
  } else if (n_mean_param == 2){
    #In this case the design matrices are an intercept and one covariate x1 which is drawn from a normal distribution with mean 0 and sd = 1.
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals_in_cluster), rnorm(n_individuals_in_cluster))), simplify = F)
    x_1 <- c(sapply(design_matrices, function(x) x[,2]))
    x_2 <- NA
    
  } else {
    #In this case the design matrices are an intercept and one covariate x1 which is drawn from a normal distribution with mean 0 and sd = 1.
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals_in_cluster), 
                                                             rnorm(n_individuals_in_cluster),
                                                             rnorm(n_individuals_in_cluster))), simplify = F)
    x_1 <- c(sapply(design_matrices, function(x) x[,2]))
    x_2 <- c(sapply(design_matrices, function(x) x[,3]))
  }
   #replicate(n_clusters, as.matrix(rep(1, n_individuals_in_cluster)), simplify = F)
  
  #---------------------------------------------------------------------------------------------------------------------------------
  
  
  
  
  #--------------------Generate data based on covariance and mean value structure and specified parameter values--------------------
  
  
  # Sampling n_clusters outcomes from a multivariate normal distribution with specified mean and covariance structure with n_individuals in each cluster
  # Adding mean-value parameters if specified
  
  if (n_cov_param == 1){
    if (n_mean_param == 1){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0)))
    } else if (n_mean_param == 2){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0))) + beta_1 * x_1
    } else {
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0))) + beta_1 * x_1 + beta_2 * x_2
    }
  } else if (n_cov_param == 2){
    if (n_mean_param == 1){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1)))
    } else if (n_mean_param == 2){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1))) + beta_1 * x_1
    } else {
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1))) + beta_1 * x_1 + beta_2 * x_2
    }
  } else if (n_cov_param == 3) {
    if (n_mean_param == 1){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2)))
    } else if (n_mean_param == 2){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2))) + beta_1 * x_1
    } else {
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2))) + beta_1 * x_1 + beta_2 * x_2
    }
  } else {
    if (n_mean_param == 1){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2 + gamma3_matrix * sigma_3)))
    } else if (n_mean_param == 2){
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2 + gamma3_matrix * sigma_3))) + beta_1 * x_1
    } else {
      y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2 + gamma3_matrix * sigma_3))) + beta_1 * x_1 + beta_2 * x_2
    }
  }
  
 # y <- c(replicate(n_clusters, mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), 
 #                                 Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2))) + beta_1 * x_1 + beta_2 * x_2
  
  # Collecting simulated data in dataframe
  DF <- data.frame(y = y, klasse = klasser, subklasse = subklasser, anden = other, x1 = x_1, x2 = x_2)
  
  # Collecting outcomes in a list of length n_clusters containing lists of outcomes each of length n_individuals_in_cluster
  outcome_list <- split(DF$y, rep(1:n_clusters, each = n_individuals_in_cluster, length.out = length(DF$y)))
  
  
  
  return(list('design_matrices' = design_matrices, 'semi_def_matrices' = semi_def_matrices, 'outcome_list' = outcome_list, 'DF' = DF))
}




