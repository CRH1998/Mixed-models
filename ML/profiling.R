#Profiling algorithm

#Define number of clusters and number of individuals in each cluster
n_clusters = 1000
n_individuals_in_cluster = 10


#Generate large dataset
large_dataset <- large_dataset_generator(n_clusters = n_clusters, n_individuals_in_cluster = n_individuals_in_cluster, seed = 1)


#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list

#Run and profile ML-algorithm

Rprof()
Rprof(NULL)

find_mle_parameters(init_params = c(1,1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome)


omega_func(semi_def_matrix_list, sigma2_vec)
omega_func_test(semi_def_matrix_list, sigma2_vec)

summaryRprof()
