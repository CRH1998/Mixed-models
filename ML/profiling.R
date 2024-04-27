#Profiling algorithm

#Define number of clusters and number of individuals in each cluster
n_clusters = 100000
n_individuals_in_cluster = 10


#Generate large dataset
large_dataset <- large_dataset_generator(n_clusters = n_clusters, 
                                         n_individuals_in_cluster = n_individuals_in_cluster, 
                                         sigma_0 = 1,
                                         sigma_1 = 3,
                                         seed = 1)


#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list

DF <- large_dataset$DF

model <- lme4::lmer(y ~ 1 + (1|subklasse), data=DF, REML = T)
summary_model <- summary(model)

#Run and profile ML-algorithm

find_mle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 0.1, tolerance = 1e-6)

system.time(find_remle_parameters(init_params = c(1,1), design_matrices = design_matrices, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 0.1, tolerance = 1e-6))
find_remle_parameters(init_params = c(1,1), design_matrices = design_matrices, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-12)

#40.477   0.947  41.453 

Rprof()
summaryRprof()
Rprof(NULL)
