#Profiling algorithm

#Define number of clusters and number of individuals in each cluster
n_clusters = 1000
n_individuals_in_cluster = 20


#Generate large dataset
large_dataset <- dataset_generator(n_clusters = n_clusters, 
                                         n_individuals_in_cluster = n_individuals_in_cluster, 
                                         mean_val = 3,
                                         sigma_0 = 0.000000000005,
                                         sigma_1 = 0.5,
                                         sigma_2 = 0.1,
                                         seed = NA)


#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list

DF <- large_dataset$DF
#View(DF)


#Run and profile ML-algorithm

find_mle_parameters(init_params = c(1,1,1,1), design_matrices = design_matrices, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-12)

find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-6, small_value_threshold = 1e-12)

system.time(find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-12))

#bruger   system forl?bet 
#1.81     0.00     1.81 


start <- Sys.time() 
find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-12, small_value_threshold = 0.1)
end <- Sys.time()

end - start

Rprof()
find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-12)
summaryRprof()
Rprof(NULL)




model <- lme4::lmer(y ~ 1 + (1|klasse) + (1|subklasse), data=DF, REML = T)
summary_model <- summary(model)
summary_model












