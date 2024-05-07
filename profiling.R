#Profiling algorithm

#Define number of clusters and number of individuals in each cluster
n_clusters = 400000
n_individuals_in_cluster = 16


#Generate large dataset
large_dataset <- dataset_generator(n_clusters = n_clusters, 
                                         n_individuals_in_cluster = n_individuals_in_cluster, 
                                         mean_val = 1,
                                         beta_1 = 5,
                                         beta_2 = 10,  
                                         sigma_0 = 3,
                                         sigma_1 = 5,
                                         sigma_2 = 7,
                                         seed = NA)


#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list

DF <- large_dataset$DF
#View(DF)


#Run and profile ML-algorithm

find_mle_parameters(init_params = c(2,2,2,2,2,2), design_matrices = design_matrices, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1)

find_remle_parameters(init_params = c(2,2,2), design_matrices = design_matrices, 
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1)

system.time(find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1))





start <- Sys.time() 
find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1)
end <- Sys.time()

end - start

#RMLE run time for 400000 clusters and 16 individuals in each clusters and 3 variance components and three fixed components: 9,4 mins


Rprof()
find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, 
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-12)
summaryRprof()
Rprof(NULL)




model <- lme4::lmer(y ~ 1 + (1|klasse) + (1|subklasse), data=DF, REML = T)
summary_model <- summary(model)
summary_model












