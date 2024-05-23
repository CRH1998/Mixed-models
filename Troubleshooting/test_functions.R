

#Define number of clusters and number of individuals in each cluster
n_clusters = 10000
n_individuals_in_cluster = 10


#Generate large dataset
large_dataset <- large_dataset_generator(n_clusters = n_clusters, sigma_0 = 3, sigma_1 = 3,
                                         n_individuals_in_cluster = n_individuals_in_cluster, seed = 1)

#large_dataset$semi_def_matrices[[1]]
#Extract dataframe from generated data
data <- large_dataset$DF
#View(data)

#Run and time lmer() function
#start <- Sys.time()
model <- lme4::lmer(y ~ 1 + (1|subklasse), data=data, REML=T)
summary_model <- summary(model)
summary_model
#end <- Sys.time()

#paste0('R function time: ', end - start)

#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list

#Run and time ML-algorithm
start <- Sys.time()
find_mle_parameters(init_params = c(1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome)
end <- Sys.time()

paste0('Own function time: ', end - start)


#4.35145497322083
#3.97797608375549

































mat1 <- matrix(data = rep(1, 100), nrow = 10)
mat2 <- matrix(data = rep(2, 100), nrow = 10)
mat3 <- matrix(data = rep(3, 100), nrow = 10)
mat4 <- matrix(data = rep(4, 100), nrow = 10)
mat5 <- matrix(data = rep(5, 100), nrow = 10)
mat6 <- matrix(data = rep(6, 100), nrow = 10)
mat7 <- matrix(data = rep(7, 100), nrow = 10)
mat8 <- matrix(data = rep(8, 100), nrow = 10)
mat9 <- matrix(data = rep(9, 100), nrow = 10)
mat10 <- matrix(data = rep(10, 100), nrow = 10)


mat_list <- list(mat1, mat2, mat3, mat4,mat5,mat6,mat7,mat8,mat9,mat10)


#matrix(unlist(lapply(mat_list, FUN = tr_multiply_list_by_matrix, list = mat_list)), nrow = 27)


loop_function <- function(mat_list){
  
  S <- matrix(data = NA, nrow = length(mat_list), ncol = length(mat_list))
  
  for (i in 1:length(mat_list)){
    for (j in 1:length(mat_list)){
      S[i,j] <- sum(mat_list[[i]] * mat_list[[j]])
    }
  }
  return(S)
}

lapply_function <- function(mat_list){
  num_matrices <- length(mat_list)
  result <- matrix(NA, nrow = num_matrices, ncol = num_matrices)
  
  # Initialize cluster
  cl <- makeCluster(detectCores())
  
  # Splitting the task and performing parallel computation
  tasks <- vector("list", length = num_matrices)
  for (i in 1:num_matrices) {
    tasks[[i]] <- clusterApply(cl, mat_list, function(x) sapply(x, function(y) sum(mat_list[[i]] * y)))
  }
  
  # Collecting and organizing results
  for (i in 1:num_matrices) {
    result[i, ] <- unlist(unname(clusterEvalQ(cl, tasks[[i]])))
  }
  
  # Stop cluster
  stopCluster(cl)
  
  return(result)
}

lapply_function(mat_list)
loop_function(mat_list)

microbenchmark(loop_function(mat_list),lapply_function(mat_list))
microbenchmark(for (i in 1:length(mat_list)){sum(mat1 * mat_list[[i]])}, unlist(tr_multiply_list_by_matrix(mat1, mat_list)))










#----------------------------------------------------------------------------
#           Checking REML and ML codes
#----------------------------------------------------------------------------

n_clusters = 10
n_individuals_in_cluster = 10

#Generate large dataset
family_dataset <- family_dataset_generator(n_clusters = n_clusters, n_individuals_in_cluster = n_individuals_in_cluster,
                                           n_mean_param = 2, mean_val = 1, beta_1 = 5,
                                           sigma_0 = 2, sigma_1 = 1, sigma_2 = 1,
                                           seed = NA)


#Extract relevant lists from data to run FS-algorithms
design_matrices <- family_dataset$design_matrices
semi_def_matrices <- family_dataset$semi_def_matrices
outcome_list <- family_dataset$outcome_list


params <- c(1,1,1)
X <- design_matrices[[1]]
V_i <- semi_def_matrices[[1]]
y <- outcome[[1]]

# Calculating omega inverse
omega <- omega_func(semi_def_matrix = V_i, sigma2_vec = params)

# Setting very small values to 0
omega[omega < small_value_threshold] <- 0

# Adding small value to diagonal if diagonal values are very small
omega <- omega + (diag(omega) < small_value_threshold) * add_small_constant * diag(length(diag(omega)))

# Inverting omega
omega_inv <- chol2inv(chol(omega))

# Calculating P matrix
V_inv <- omega_inv

P_test <- omega_inv - omega_inv %*% X %*% solve((t(X) %*% omega_inv %*% X)) %*% t(X) %*% omega_inv
P <- P_func(omega_inv = V_inv, design_matrix = X)

P_test - P

#Multiplying the P matrix with each semi-definite matrix (variance component matrix) to save computations
A <- multiply_list_by_matrix(P, V_i)


#-------------Calculating S matrix-----------------
S <- S_matrix_reml_function(P = P, semi_def_matrix = V_i)

S

0.5 * tr(A[[1]] %*% A[[2]])

0.5 * tr(P%*%P)
0.5 * tr(P %*% V_i[[1]] %*% P %*% V_i[[2]])
0.5 * tr(P %*% V_i[[1]] %*% P %*% V_i[[3]])
0.5 * tr(P %*% V_i[[2]] %*% P %*% V_i[[1]])
0.5 * tr(P %*% V_i[[2]] %*% P %*% V_i[[2]])
0.5 * tr(P %*% V_i[[2]] %*% P %*% V_i[[3]])
0.5 * tr(P %*% V_i[[3]] %*% P %*% V_i[[1]])
0.5 * tr(P %*% V_i[[3]] %*% P %*% V_i[[2]])
0.5 * tr(P %*% V_i[[3]] %*% P %*% V_i[[3]])


#-------------Calculating scores-------------------
score <- reml_score_func(P = P, outcomes = y, semi_def_matrix = V_i)


0.5*(-tr(P) + t(y) %*% P %*% P %*% y)
0.5*(-tr(P %*% V_i[[2]]) + t(y) %*% P %*% V_i[[2]] %*% P %*% y)
0.5*(-tr(P %*% V_i[[3]]) + t(y) %*% P %*% V_i[[3]] %*% P %*% y)



init_params <- c(1,1,1)
small_value_threshold = 1e-12
add_small_constant = 1e-12

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











