

#Define number of clusters and number of individuals in each cluster
n_clusters = 100
n_individuals_in_cluster = 10


#Generate large dataset
large_dataset <- large_dataset_generator(n_clusters = n_clusters, n_individuals_in_cluster = n_individuals_in_cluster, seed = 1)


#Extract dataframe from generated data
data <- large_dataset$DF


#Run and time lmer() function
start <- Sys.time()
model <- lme4::lmer(y ~ 1 + (1|klasse) + (1|subklasse), data=data, REML=F)
summary_model <- summary(model)
summary_model
end <- Sys.time()

paste0('R function time: ', end - start)

#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list

#Run and time ML-algorithm
start <- Sys.time()
find_mle_parameters(init_params = c(1,1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome)
end <- Sys.time()

paste0('Own function time: ', end - start)



