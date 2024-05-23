######################################################
#                                                    #
# This code constructs a test family dataset for     #
# testing the fisher scoring algorithms              #
#                                                    #
######################################################

#Construct function that pairs parents with children
pair_parents_with_children <- function(parents, unassigned_individuals){
  
  sub_family_df <- data.frame()
  
  for (i in 1:nrow(parents)){
    # Assign potential children to parents
    no_children <- rbinom(1,nrow(unassigned_individuals),0.25)
    
    children <- unassigned_individuals[sample(nrow(unassigned_individuals), no_children),]

    
    # Construct sub_family dataframe
    ifelse(no_children == 0, #If the number of sampled children is 0, the sub_family dataframe is empty
           sub_family_df <- rbind(sub_family_df, data.frame(matrix(ncol = 4, nrow = 0))), 
           sub_family_df <- rbind(sub_family_df, data.frame(children, parents[i,1], parents[i,3]))
           )
    
    #Remove children who were assigned in previous statement
    unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% children$id),]
  }
  
  names(sub_family_df) <- c("id", "sex","dadid", "momid")
  return(list('sub_family_df' = sub_family_df, 'unassigned_individuals' = unassigned_individuals))
}


# Construct function that pairs children with spouses
pair_children_with_spouses <- function(family_df, unassigned_individuals){
  
  # Get sex of children
  male_children <- family_df[family_df$sex == 1,]
  female_children <- family_df[family_df$sex == 2,]
  
  # Get ids of men and women in the id column
  unassigned_men <- unassigned_individuals[unassigned_individuals$sex == 1,]
  unassigned_women <- unassigned_individuals[unassigned_individuals$sex == 2,]
  
  # Sample spouses for children among the remaining
  no_male_spouses <- sample(0:min(nrow(female_children), nrow(unassigned_men)), 1)
  no_female_spouses <- sample(0:min(nrow(male_children), nrow(unassigned_women)), 1)
  
  married_male_unassigned <- unassigned_men[sample(nrow(unassigned_men), no_male_spouses),]
  married_female_unassigned <- unassigned_women[sample(nrow(unassigned_women), no_female_spouses),]
  
  # Sample children:
  married_male_children <- male_children[sample(nrow(male_children), no_female_spouses),][,c("id", "sex")]
  married_female_children <- female_children[sample(nrow(female_children), no_male_spouses),][,c("id", "sex")]
  
  # Construct parent dataframe
  parents <- data.frame(rbind(cbind(married_male_children, married_female_unassigned), 
                              cbind(married_male_unassigned, married_female_children)))
  
  # Update assigned individuals
  ifelse(no_male_spouses == 0 & no_female_spouses == 0, 
         assigned_individuals <- data.frame(matrix(ncol = 4, nrow = 0)),
         assigned_individuals <- data.frame(rbind(married_male_unassigned, married_female_unassigned), 0, 0))
  
  names(assigned_individuals) <- c("id", "sex","dadid", "momid")
  
  # Update unassigned individuals
  unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% assigned_individuals$id),]
  
  return(list('parents' = parents, 'unassigned_individuals' = unassigned_individuals, 'assigned_individuals' = assigned_individuals))
}

# Construct function that generates kinship matrix
generate_kinship_matrix <- function(n_individuals_in_cluster, seed = NA){
  
  if (!is.na(seed)){
    set.seed(seed)
  }
  

  # Specify id of each individual in the family
  id  <- 1:n_individuals_in_cluster
  
  # Randomly assign sex to each individual in the family
  sex <- sample(c(1,2), n_individuals_in_cluster, prob = c(1/2,1/2), replace = T)
  
  # Unassigned individuals
  unassigned_individuals <- data.frame(id,sex)
  
  # Get ids of men and women in the id column
  men <- unassigned_individuals[unassigned_individuals$sex == 1,]
  women <- unassigned_individuals[unassigned_individuals$sex == 2,]
  
  # Pair random man with random woman to create first parents
  first_dad <- men[sample(nrow(men), 1),]
  first_mom <- women[sample(nrow(women), 1),]
  
  parents <- data.frame(first_dad, first_mom)
  
  # Remove first parents from unassigned individuals
  unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% c(first_dad$id, first_mom$id)),]
  
  full_family_df <- data.frame(rbind(first_dad, first_mom), 0, 0)
  names(full_family_df) <- c("id", "sex","dadid", "momid")
  
  while(nrow(unassigned_individuals) > 0){
    
    sub_family <- pair_parents_with_children(parents, unassigned_individuals)
      
    # Define sub_family_df
    sub_family_df <- sub_family$sub_family_df
      
    # Update unassigned individuals
    unassigned_individuals <- sub_family$unassigned_individuals

    # Assign spouses to children
    parents_list <- pair_children_with_spouses(sub_family_df, unassigned_individuals)
      
    # Update potential parents
    if (nrow(parents_list$parents) > 0){
      parents <- parents_list$parents
    }
    
    # Update assigned individuals, that is unassigned individuals who were married to one of the children
    assigned_individuals <- parents_list$assigned_individuals
    unassigned_individuals <- parents_list$unassigned_individuals
    full_family_df <- rbind(full_family_df, sub_family_df, assigned_individuals)
    
  }
  
  full_family_df <- full_family_df %>% arrange(id)
  
  kinship_matrix <- kinship(id = full_family_df$id, 
                            dadid = full_family_df$dadid, 
                            momid = full_family_df$momid, 
                            sex = full_family_df$sex)
  
  return(kinship_matrix)
}

# Construct function that generates multiple kinship matrices
generate_n_kinship_matrices <- function(n_individuals_in_cluster, n_clusters, seed = NA){
  kinship_matrices <- list()
  
  for (i in 1:n_clusters){
    kinship_matrices[[i]] <- generate_kinship_matrix(n_individuals = n_individuals_in_cluster)
  }
  
  return(kinship_matrices)
}








#-----------------------------------------------
#   Generate family data for testing functions
#-----------------------------------------------


# Generate data set for testing
family_dataset_generator <- function(n_clusters, n_individuals_in_cluster, n_mean_param = 3,
                                     mean_val = 1, beta_1 = 3, beta_2 = 3, sigma_0 = 1, sigma_1 = 1, sigma_2 = 1, seed = 1){
  
  # Set seed if specified for sampling later on
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  
  
  #-------------------------------------------Generate covariance structure-------------------------------------------
  
  #Diagonal matrix with ones in the diagonal and dimension corresponding to the number of individuals in each cluster (intercept)
  gamma0_matrix <- as.matrix(diag(1, nrow = n_individuals_in_cluster))
  
  
  #Matrix of ones indicating family relation
  gamma1_matrix <- matrix(1, nrow = n_individuals_in_cluster, ncol = n_individuals_in_cluster) 
  #gamma1_matrix <- as.matrix(bdiag(matrix(1, nrow = floor(n_individuals_in_cluster/2), ncol = floor(n_individuals_in_cluster/2)), matrix(1, nrow = ceiling(n_individuals_in_cluster/2), ncol = ceiling(n_individuals_in_cluster/2))))
  

  
  
  # Defining semi_def_matrices which is a list containing n_clusters replicas of the gamma_list
  semi_def_matrices <- replicate(n_clusters,
                                 list(gamma0_matrix, gamma1_matrix, generate_kinship_matrix(n_individuals = n_individuals_in_cluster)*2),
                                 simplify = F)
  
  #---------------------------------------------------------------------------------------------------------------------------------
  
  
  
  
  #-------------------------------------------Construct mean value structure---------------------------------------------------------
  # Defining the design_matrices 
  if (n_mean_param == 1){
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals_in_cluster))), simplify = F)
  } else if (n_mean_param == 2){
    #In this case the design matrices are an intercept and one covariate x1 which is drawn from a normal distribution with mean 0 and sd = 1.
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals_in_cluster), rnorm(n_individuals_in_cluster))), simplify = F)
  } else {
    #In this case the design matrices are an intercept and one covariate x1 which is drawn from a normal distribution with mean 0 and sd = 1.
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals_in_cluster), 
                                                             rnorm(n_individuals_in_cluster),
                                                             rnorm(n_individuals_in_cluster),
                                                             rnorm(n_individuals_in_cluster),
                                                             rnorm(n_individuals_in_cluster),
                                                             rnorm(n_individuals_in_cluster))), simplify = F)
  }
  #---------------------------------------------------------------------------------------------------------------------------------
  
  
  
  
  #--------------------Generate data based on covariance and mean value structure and specified parameter values--------------------
  
  
  # Sampling n_clusters outcomes from a multivariate normal distribution with specified mean and covariance structure with n_individuals in each cluster
  # Adding mean-value parameters if specified
  outcome_list <- list()
  if (n_mean_param == 1){
    for (i in 1:n_clusters){
      outcome_list[[i]] <- mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = semi_def_matrices[[i]][[1]] * sigma_0 + semi_def_matrices[[i]][[2]] * sigma_1 + semi_def_matrices[[i]][[3]] * sigma_2)
    }
  } else if (n_mean_param == 2){
    for (i in 1:n_clusters){
      outcome_list[[i]] <- mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = semi_def_matrices[[i]][[1]] * sigma_0 + semi_def_matrices[[i]][[2]] * sigma_1 + semi_def_matrices[[i]][[3]] * sigma_2) + beta_1 * design_matrices[[i]][,2]
    }
  } else {
    for (i in 1:n_clusters){
      outcome_list[[i]] <- mvrnorm(n = 1, mu = rep(mean_val, n_individuals_in_cluster), Sigma = semi_def_matrices[[i]][[1]] * sigma_0 + semi_def_matrices[[i]][[2]] * sigma_1 + semi_def_matrices[[i]][[3]] * sigma_2)  
      + beta_1 * design_matrices[[i]][,2] + beta_2 * design_matrices[[i]][,3] + 4 * design_matrices[[i]][,4] + 5 * design_matrices[[i]][,5] + 6 * design_matrices[[i]][,6]
    }
  }
  return(list('design_matrices' = design_matrices, 'semi_def_matrices' = semi_def_matrices, 'outcome_list' = outcome_list))
}











