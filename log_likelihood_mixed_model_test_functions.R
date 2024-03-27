


func_input <- gen_func_input(n_blokke = 3)

design_matrices <- func_input$design_matrices
semi_def_matrices <- func_input$semi_def_matrices
semi_def_matrix <- semi_def_matrices[[1]]
outcomes <- func_input$outcome_list
summary <- func_input$summary_model
summary


sigma2_vec <- c(2.71, 17.73)
mean_value_param <- 8.321
parameters <- c(mean_value_param, sigma2_vec)

omega_matrix <- omega_func(semi_def_matrix, sigma2_vec)
omega_inv <- chol2inv(chol(omega_matrix))




find_mle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcomes)

find_remle_parameters(init_params = c(1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcomes)


microbenchmark(find_mle_parameters(init_params = c(1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcomes))



list_outer <- function(list1, list2){
  mapply(element_product, list1, list2)
}



element_product2 <- function(elem1, elem2){
  return(0.5 * tr(elem1 %*% elem2))
}

element_product3 <- function(elem1, elem2){
  return(elem1 %*% elem2)
}




S_elements <- function(matrix1, matrix2){
  return(0.5 * tr(matrix1 %*% matrix2))
}

S <- function(a,b, element_product1) {
  outer(a, b, function(x,y) vapply(seq_along(x), function(i) element_product(x[[i]], y[[i]]), numeric(1)))
}



f2 <- function(semi_def_matrix, fixed_matrix, element_product2){
  A <- multiply_list_by_matrix(fixed_matrix, semi_def_matrix)
  
  outer(A, A, function(x,y) vapply(seq_along(x), function(i) element_product2(x[[i]], y[[i]]), numeric(1)))
}


my_array <- array(unlist(semi_def_matrix), dim = c(nrow(semi_def_matrix[[1]]), ncol(semi_def_matrix[[1]]), length(semi_def_matrix)))
outer(my_array, my_array, FUN = element_product3)


f1(semi_def_matrix, semi_def_matrix, omega_inv, element_product)
f2(semi_def_matrix, omega_inv, element_product2)

microbenchmark(f1(semi_def_matrix, semi_def_matrix, omega_inv, element_product), f2(semi_def_matrix, omega_inv, element_product2))


microbenchmark(f1(semi_def_matrix, semi_def_matrix, omega_inv, element_product), 
               S_matrix_function(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv))

microbenchmark(S_matrix_function_test(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv), 
               S_matrix_function(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv))


microbenchmark(multiply_list_by_matrix(omega_inv, semi_def_matrix))





























##########################################
#               REML score
##########################################

#-------------Generating function input---------------
func_input <- gen_func_input(n_blokke = 3, REML = T)


#-------------Design matrices-------------
design_matrices <- func_input$design_matrices
design_matrices


#-------------Semi definite matrices (tau matrices)-------------
semi_def_matrices <- func_input$semi_def_matrices
semi_def_matrices


#-------------List of outcomes-------------
outcome_list <- func_input$outcome_list
outcome_list


#-------------Summary of model------------------------
summary_model <- func_input$summary_model
summary_model



#-------------Initial parameters-------------
sigma2_vec <- c(20,2) #True parameters: c(26.925, 2.709)





######################################################
#         Function inputs for single one block
######################################################

X <- design_matrices[[1]]
semi_def_matrix <- semi_def_matrices[[1]]
y <- outcome_list[[1]]


#---------Variance components----------------
I_matrix <- semi_def_matrix[[1]]
A_matrix <- semi_def_matrix[[2]]


#--------Covariance matrix-------------------
omega <- omega_func(semi_def_matrix = semi_def_matrix, sigma2_vec = sigma2_vec)
omega_inv <- chol2inv(chol(omega))



#---------------P matrix---------------------
P <- P_func(omega_inv, X)


#---------------Scores---------------
scores <- reml_score_func(P = P, semi_def_matrix = semi_def_matrix, outcomes = y)


- 0.5 * tr(P %*% I_matrix) + 0.5 * t(y) %*% P %*% I_matrix %*% P %*% y
- 0.5 * tr(P %*% A_matrix) + 0.5 * t(y) %*% P %*% A_matrix %*% P %*% y


#---------------S matrix---------------
S_matrix <- S_matrix_reml_function(semi_def_matrix = semi_def_matrix, P = P)
S_matrix_ML <- S_matrix_function(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv)


0.5 * tr(P %*% I_matrix %*% P %*% I_matrix)
0.5 * tr(P %*% I_matrix %*% P %*% A_matrix)
0.5 * tr(P %*% A_matrix %*% P %*% I_matrix)
0.5 * tr(P %*% A_matrix %*% P %*% A_matrix)



sigma2_vec + solve(S_matrix) %*% scores



reml_score_fisher_function(X, semi_def_matrix = semi_def_matrix, outcomes = y, params = sigma2_vec)
find_remle_parameters(sigma2_vec, design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome_list)







P <- P_func(omega_inv, X)
Py <- Py_func_P(P, y)
y_t_P <- y_t_P_func(P, y)

t(S_matrix_reml_function(semi_def_matrix, P))


P %*% semi_def_matrix[[2]]
semi_def_matrix[[2]]

omega_inv - omega_inv %*% X %*% solve(t(X) %*% omega_inv %*% X) %*% t(X) %*% omega_inv




solve(S_matrix_reml_function(semi_def_matrix, P))

model_summary


P %*% outcome
Py_func(omega_inv = omega_inv, design_matrix = X, outcomes = outcome, beta)
















####################################
#       Generate large dataset     #
####################################

n_clusters = 10
n_individuals_in_cluster = 10






large_dataset <- large_dataset_generator(5000, 10, seed = 1)



data <- large_dataset$DF

start <- Sys.time()
model <- lme4::lmer(y ~ 1 + (1|klasse) + (1|subklasse), data=data, REML=F)
summary_model <- summary(model)
summary_model
end <- Sys.time()

paste0('R function time: ', end - start)

design_matrices <- large_dataset$design_matrices
semi_def_matrices <- large_dataset$semi_def_matrices
outcome <- large_dataset$outcome_list



start <- Sys.time()
find_mle_parameters(init_params = c(1,1,1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome)
end <- Sys.time()

paste0('Own function time: ', end - start)


#Measuring time

#n_clusters:          10              100               500             1000            10000               50000              100         
#n_individuals:       10              100               100             100             10                  20                 750
#time:              0.401541 sec    0.849829 sec      3.565483 sec    6.822527 sec    7.28249 sec       41.53817 secs     3.870443 mins

find_remle_parameters(init_params = c(1,1), design_matrices = design_matrices, semi_def_matrices = semi_def_matrices, outcome_list = outcome)








#Optimize omega matrix calculation
n_clusters = 1000
n_individuals_in_cluster = 300
large_dataset <- large_dataset_generator(n_clusters, n_individuals_in_cluster, seed = 1)

sigma2_vec <- c(1.5,2,2.5)
semi_def_matrices <- large_dataset$semi_def_matrices

semi_def_matrix_list <- semi_def_matrices[[1]]




microbenchmark(omega_func(semi_def_matrix_list, sigma2_vec), omega_func_test(semi_def_matrix_list, sigma2_vec))









