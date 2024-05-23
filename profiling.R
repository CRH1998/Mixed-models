#Profiling algorithm

#Define number of clusters and number of individuals in each cluster
n_clusters = 400
n_individuals_in_cluster = 12


#Generate large dataset
large_dataset <- dataset_generator(n_clusters = n_clusters, n_cov_param = 2, n_mean_param = 1,
                                         n_individuals_in_cluster = n_individuals_in_cluster, 
                                         mean_val = 1,
                                         beta_1 = 5,
                                         beta_2 = 3,  
                                         sigma_0 = 3,
                                         sigma_1 = 5,
                                         sigma_2 = 7,
                                         sigma_3 = 3,
                                         seed = NA)


#Extract relevant lists from data to run ML-algorithm
design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list



#Run and profile ML-algorithm
mle_parameters <- find_mle_parameters(init_params = c(2,2,2), design_matrices = design_matrices, tolerance = 1e-9, 
                    semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1)


rmle_parameters <- find_remle_parameters(init_params = c(2,2), design_matrices = design_matrices, tolerance = 1e-12,
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1)



#Timing REML FS-algorithm
start <- Sys.time()
find_remle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1)
end <- Sys.time()
end - start

#RMLE run time for 400000 clusters and 16 individuals in each clusters and 3 variance components and 3 fixed components robust approach: 9,4 mins
#RMLE run time for 400000 clusters and 16 individuals in each clusters and 3 variance components and 3 fixed components semi-robust approach: 7,9 mins
#RMLE run time for 400000 clusters and 16 individuals in each clusters and 3 variance components and 3 fixed components non-robust approach: 7,6 mins




#Profiling REML FS-algorithm
Rprof()
find_remle_parameters(init_params = c(1,1,1,1), design_matrices = design_matrices, 
                      semi_def_matrices = semi_def_matrices, outcome_list = outcome, update_step_size = 1, tolerance = 1e-12)
summaryRprof()
Rprof(NULL)


#Testing log-likelihood functions
log_likelihood(design_matrices, semi_def_matrices, outcome, c(0.93232, 4.99647, 3.03114, 2.831, 4.130, 7.545, 2.963))


restricted_log_likelihood(design_matrices, semi_def_matrices, outcome, c(3.090583,4.671067,6.638998))

lower_bounds <- c(-Inf, -Inf, -Inf, 0.000001, 0.000001, 0.000001, 0.000001)
upper_bounds <- c(Inf, Inf, Inf, Inf, Inf, Inf, Inf)

optim(par = c(1,1,1,1,1,1,1), 
      fn = log_likelihood, 
      design_matrices = design_matrices, 
      semi_def_matrices = semi_def_matrices, 
      outcome_list = outcome,
      lower = lower_bounds,
      upper = upper_bounds,
      method = 'L-BFGS-B',
      control=list(fnscale=-1))

#Testing lme4::lmer:
DF <- large_dataset$DF
#View(DF)

model <- lme4::lmer(y ~ 1 + (1|klasse), data=DF, REML = T)
summary_model <- summary(model)
summary_model$logLik
log_likelihood(design_matrices, semi_def_matrices, outcome, mle_parameters)
summary_model












