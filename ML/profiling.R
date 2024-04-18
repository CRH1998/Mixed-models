#Profiling algorithm

#Define number of clusters and number of individuals in each cluster
n_clusters = 100000
n_individuals_in_cluster = 50


#Generate large dataset
large_dataset <- large_dataset_generator(n_clusters = n_clusters, n_individuals_in_cluster = n_individuals_in_cluster, seed = 1)


#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list

DF <- large_dataset$DF
model <- lme4::lmer(y ~ 1 + (1|klasse) + (1|subklasse), data=DF, REML=F)
summary_model <- summary(model)

#Run and profile ML-algorithm

find_mle_parameters(init_params = c(1,1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1)



#find_remle_parameters(init_params = c(1,1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 0.1)