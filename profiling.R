#Profiling algorithm

#Define number of clusters and number of individuals in each cluster
n_clusters = 1
n_individuals_in_cluster = 24


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


restricted_log_likelihood_block(design_matrices[[1]], semi_def_matrices[[1]], outcome[[1]], c(2.405801, 5))

lower_bounds <- c(0.000001, 0.000001)
upper_bounds <- c(Inf, Inf)

optim(par = c(1,2), 
      fn = restricted_log_likelihood_block,
      design_matrix = design_matrices[[1]], 
      semi_def_matrix = semi_def_matrices[[1]], 
      outcomes = outcome[[1]],
      lower = lower_bounds,
      upper = upper_bounds,
      method = 'L-BFGS-B',
      control=list(fnscale=-1))









#Define number of clusters and number of individuals in each cluster
n_clusters = 1
n_individuals_in_cluster = 96


#Generate large dataset
large_dataset <- dataset_generator(n_clusters = n_clusters, n_cov_param = 2, n_mean_param = 1,
                                   n_individuals_in_cluster = n_individuals_in_cluster, 
                                   mean_val = 0,
                                   beta_1 = 0,
                                   beta_2 = 0,  
                                   sigma_0 = 10,
                                   sigma_1 = 10,
                                   sigma_2 = 7,
                                   sigma_3 = 3,
                                   seed = NA)

X <- large_dataset$design_matrices
semi_def_matrix <- large_dataset$semi_def_matrices
y <- large_dataset$outcome_list

get_outcome_variance(y)
get_outcome_mean(y)

#Testing lme4::lmer:
DF <- large_dataset$DF


model <- lme4::lmer(y ~ 1 + (1|klasse), data=DF, REML = T)
summary_model <- summary(model)

lme4:::devCrit
lme4:::devcomp
model@devcomp

summary_model
summary_model$logLik
log_likelihood(design_matrices = X,semi_def_matrices = semi_def_matrix,outcome_list = y,parameters = c(-3.352,10.785,3.953))
restricted_log_likelihood(design_matrices = X,semi_def_matrices = semi_def_matrix,outcome_list = y,parameters = c(-3.352,10.785, 8.131))



lower_bounds <- c(0.000001, 0.000001)
upper_bounds <- c(Inf, Inf)

optim(par = c(1,2), 
      fn = restricted_log_likelihood,
      design_matrices = X, 
      semi_def_matrices = semi_def_matrix, 
      outcome_list = y,
      lower = lower_bounds,
      upper = upper_bounds,
      method = 'L-BFGS-B',
      control=list(fnscale=-1))




n_clusters = 2
n_individuals_in_cluster = 100

familiy_dataset <- family_dataset_generator(n_clusters, n_individuals_in_cluster, 
                                            n_variance_param = 2, n_mean_param = 1,
                                            mean_val = 0, sigma_0 = 4, sigma_1 = 4, seed = 1)

X <- familiy_dataset$design_matrices
semi_def_matrix <- familiy_dataset$semi_def_matrices
y <- familiy_dataset$outcome_list

get_outcome_variance(y)
get_outcome_mean(y)

restricted_log_likelihood(design_matrices = X,semi_def_matrices = semi_def_matrix,outcome_list = y,parameters = c(1,1))

lower_bounds <- c(0.000001, 0.000001)
upper_bounds <- c(Inf, Inf)

optim(par = c(1,2), 
      fn = restricted_log_likelihood,
      design_matrices = familiy_dataset$design_matrices, 
      semi_def_matrices = familiy_dataset$design_matrices, 
      outcome_list = familiy_dataset$outcome_list,
      lower = lower_bounds,
      upper = upper_bounds,
      method = 'L-BFGS-B',
      control=list(fnscale=-1))












params <- c(3.1447^2,1.8563^2,1.6751^2)

#Extract relevant lists from data to run REML_ll
X1 <- large_dataset$design_matrices[[1]]
semi_def_matrix1 <- large_dataset$semi_def_matrices[[1]]
y1 <- large_dataset$outcome_list[[1]]

X2 <- large_dataset$design_matrices[[2]]
semi_def_matrix2 <- large_dataset$semi_def_matrices[[2]]
y2 <- large_dataset$outcome_list[[2]]


REML_p_ll_block1 <- REML_p_ll_blocks(X = X1, semi_def_matrix = semi_def_matrix1, y = y1, params = params)
REML_p_ll_block2 <- REML_p_ll_blocks(X = X2, semi_def_matrix = semi_def_matrix2, y = y2, params = params)



ytVinvy1 <- REML_p_ll_block1$ytVinvy[[1]]
XtVinvX1 <- REML_p_ll_block1$XtVinvX[[1]]
log_det_V1 <- REML_p_ll_block1$log_det_V

ytVinvy2 <- REML_p_ll_block2$ytVinvy[[1]]
XtVinvX2 <- REML_p_ll_block2$XtVinvX[[1]]
log_det_V2 <- REML_p_ll_block2$log_det_V

-0.5*((2*96-1)*log(ytVinvy1 + ytVinvy2) + log(XtVinvX1 + XtVinvX2) + log_det_V1 + log_det_V2)


REML_ll(X = X, semi_def_matrix = semi_def_matrix, y = y, params = c(2.273, 3.030, 3))
Reduce('+',Map(REML_ll, X, semi_def_matrix, y, MoreArgs = list(c(2.273, 3.030))))









