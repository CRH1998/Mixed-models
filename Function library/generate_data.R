######################################################
#                                                    #
# This code construct a test dataset for testing     #
# the fisher scoring algorithms                      #
######################################################

#-------------------------------------------
#   Generate data for testing functions
#-------------------------------------------


# Generate data set for testing
dataset_generator <- function(n_clusters, n_individuals_in_cluster, mean_val = 1, sigma_0 = 1, sigma_1 = 1, sigma_2 = 1, seed = 1){
  
  
  #Diagonal matrix with ones in the diagonal and dimension corresponding to the number of individuals in each cluster (intercept)
  gamma0_matrix <- as.matrix(diag(1, nrow = n_individuals_in_cluster))
  
  
  #Block diagonal matrix with two blocks of ones and dimension corresponding to the number of individuals in each cluster
  gamma1_matrix <- as.matrix(bdiag(matrix(1, nrow = floor(n_individuals_in_cluster/2), ncol = floor(n_individuals_in_cluster/2)), 
                                   matrix(1, nrow = ceiling(n_individuals_in_cluster/2), ncol = ceiling(n_individuals_in_cluster/2))))
  
    #Corresponding feature
    klasser <- rep(seq(1:(n_clusters*2)), each = n_individuals_in_cluster/2)
  
  
  
  #Block diagonal matrix with four blocks of ones and dimensions corresponding to the number of individuals in each cluster
  gamma2_matrix <- as.matrix(bdiag(matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4),
                                   matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4),
                                   matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4),
                                   matrix(1, nrow = n_individuals_in_cluster/4, ncol = n_individuals_in_cluster/4)))
  
    #Corresponding feature
    subklasser <- rep(seq(1:(n_clusters*4)), each = n_individuals_in_cluster/4) # rep(rep(seq(1:4),each = n_individuals_in_cluster/4), n_clusters)
  
  
  # Collecting gamma matrices into list
  gamma_list <- list(gamma0_matrix, gamma1_matrix, gamma2_matrix)
  
  
  # Defining semi_def_matrices which is a list containing n_clusters replicas of the gamma_list
  semi_def_matrices <- replicate(n_clusters, gamma_list, simplify = F)
  
  # Defining the design_matrices which is a list of n_clusters replicas of 1 vectors of length n_individuals_in_cluster
  design_matrices <- replicate(n_clusters, as.matrix(rep(1, n_individuals_in_cluster)), simplify = F)
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  
  # Sampling n_clusters outcomes from a multivariate normal distribution with n_individuals in cluster sample size
  y <- c(replicate(n_clusters, 
                   mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), 
                                  Sigma = gamma0_matrix * sigma_0 + gamma1_matrix * sigma_1 + gamma2_matrix * sigma_2)))
  
  # Collecting simulated data in dataframe
  DF <- data.frame(y = y, klasse = klasser, subklasse = subklasser)
  
  # Collecting outcomes in a list of length n_clusters containing lists of outcomes each of length n_individuals_in_cluster
  outcome_list <- split(DF$y, rep(1:n_clusters, each = n_individuals_in_cluster, length.out = length(DF$y)))
  
  
  
  return(list('design_matrices' = design_matrices, 'semi_def_matrices' = semi_def_matrices, 'outcome_list' = outcome_list, 'DF' = DF))
}




