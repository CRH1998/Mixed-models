

library(Matrix)               # For matrix manipulation
library(clusterGeneration)    # For positive semi-definite matrices
library(mvtnorm)              # For log-likelihood calculation
library(lme4)                 # For mixed model
library(microbenchmark)       # For testing function speed
library(psych)                # For calculating the trace


######################################################
#                                                    #
# Construct R-function to calculate log likelihood   #
# for block-multivariate normal gaussian where the   #
# covariance matrix is a linear combination of kno-  #
# wn semi-definite matrices                          #
#                                                    #
#                                                    #
######################################################



#Number of clusters
K <- 2

#Number of individuals in each cluster
n_i <- 3

#Number of fixed parameters
p <- 7

#Number of random parameters
m <- 5


#Design matrices
design_matrices <- list()
design_names    <- c()

# Loop to generate and add matrices to the list
for (i in 1:K) {
  # Creating a sample matrix
  matrix_data <- diag(1, n_i, p)
  
  # Adding the matrix to the list
  design_matrices[[i]] <- matrix_data
  
  # Constructing name vector
  design_names[i] <- paste0('X',i)
}
names(design_matrices) <- design_names




#Positive semi-definite matrices
semi_def_matrices <- list()
semi_def_names    <- c()
for (i in 1:K){
  partial_semi_def_matrices <- list()
  partial_semi_def_names    <- c()
  for (j in 1:(m+1)){
    # Creating a sample matrix
    semi_def_matrix <- genPositiveDefMat(dim = n_i)$Sigma
    
    # Adding the matrix to the list
    partial_semi_def_matrices[[j]] <- semi_def_matrix
    
    # Constructing name vector
    partial_semi_def_names[j] <- paste0('t',j)
  }
  names(partial_semi_def_matrices) <- partial_semi_def_names
  
  semi_def_matrices[[i]]  <- partial_semi_def_matrices
  semi_def_names[i]       <- paste0('tau',i)
  
}
names(semi_def_matrices) <- semi_def_names





outcome_list        <- list()
outcome_list_names  <- c()
for (i in 1:K){
  #Creating outcomes
  y <- rnorm(n_i, mean = 3, sd = 1)
  
  #Storing outcomes in list
  outcome_list[[i]] <- y
  
  # Constructing name vector
  outcome_list_names[i] <- paste0('y',i)
}
names(outcome_list) <- outcome_list_names


param_vec <- rnorm(p, mean = 3, sd = 1)
sigma2_vec <- rnorm(m+1,mean = 3, sd = 1)^2

parameters <- c(param_vec, sigma2_vec)













#############################################
#Testing function run time
#############################################
f <- function(par, x, y) { sum((y - (par[1]*x^2 + par[2]*x+par[3]))^2)}

x <- rnorm(100)
y <- rnorm(100, mean=.5*x^2 -2*x +4)

optim(c(1,-1,1), f, x=x, y=y)

#############################################
#Testing LL-function output
#############################################
DF <- data.frame(y = 1:15,
                 klasse=c(1,1,1,1,1,
                          2,2,2,2,2,2,
                          3,3,3,3))

model <- lme4::lmer(y ~ 1 + (1|klasse), data=DF, REML=FALSE)
summary_model <- summary(model)
LL_model <- logLik(model)
LL_model


#----------Gruppe 1-----------------
model_matrix1 <- matrix(matrix(rep(1,nrow(DF)))[1:5,])
gamma1_matrix1 <- as.matrix(diag(1, nrow = 5))
gamma2_matrix1 <- as.matrix(bdiag(matrix(1,5,5)))
outcome1 <- DF$y[1:5]

#----------Gruppe 2-----------------
model_matrix2 <- matrix(matrix(rep(1,nrow(DF)))[6:11,])
gamma1_matrix2 <- as.matrix(diag(1, nrow = 6))
gamma2_matrix2 <- as.matrix(bdiag(matrix(1,6,6)))
outcome2 <- DF$y[6:11]

#----------Gruppe 3-----------------
model_matrix3 <- matrix(matrix(rep(1,nrow(DF)))[12:15,])
gamma1_matrix3 <- as.matrix(diag(1, nrow = 4))
gamma2_matrix3 <- as.matrix(bdiag(matrix(1,4,4)))
outcome3 <- DF$y[12:15]


design_matrices <- list(model_matrix1, model_matrix2, model_matrix3)
semi_def_matrices <- list(list(gamma1_matrix1, gamma2_matrix1),
                          list(gamma1_matrix2, gamma2_matrix2),
                          list(gamma1_matrix3, gamma2_matrix3))
                          
outcome_list <- list(outcome1, outcome2, outcome3)

param_vec <- c(summary_model$coefficients[1])
sigma2_vec <- c(3.056, 7.002)
parameters <- c(param_vec, sigma2_vec)

log_likelihood(design_matrices, semi_def_matrices, outcome_list, parameters = parameters)



##########################################
#             ML-optimization
##########################################


lower_bounds <- c(-Inf, 0.000001, 0.000001)
upper_bounds <- c(Inf, Inf, Inf)

optim(par = c(1,1,1), 
      fn = log_likelihood, 
      design_matrices = design_matrices, 
      semi_def_matrices = semi_def_matrices, 
      outcome_list = outcome_list,
      lower = lower_bounds,
      upper = upper_bounds,
      method = 'L-BFGS-B',
      control=list(fnscale=-1))


#Fisher's scoring method for estimating parameters

#We wish to estimate mean vector parameter and variance/covariance parameters

#We need the Fisher information and the score wrt. each parameter

#Start with intial parameter value

#beta = 5.768514, sigma_0 = 3.055593, sigma_1 = 7.001967





semi_def_matrix <- semi_def_matrices[[1]]

beta_hat <- 1
sigma_0 <- 1
sigma_1 <- 1

omega <- bdiag(sigma_0 * semi_def_matrices[[1]][[1]] + sigma_1 * semi_def_matrices[[1]][[2]])
omega_inv <- as.matrix(solve(omega))


#Skriv en funktion, der hurtigt kan regne S-matricen ud


S_matrix_function_test <- function(semi_def_matrix, omega_inv){
  n <- length(semi_def_matrix)
  
  # Define a function to calculate the desired quantity for a given index
  calculate_S_ij <- function(i, semi_def_matrix, omega_inv) {
    term1 <- omega_inv %*% semi_def_matrix[[i]]
    0.5 * sum(diag(term1 %*% omega_inv %*% semi_def_matrix))
  }
  
  # Use lapply to calculate the result for each i
  S <- matrix(unlist(lapply(1:n, calculate_S_ij, semi_def_matrix = semi_def_matrix, omega_inv = omega_inv)), nrow = n, ncol = n)
}

S_matrix_function_test(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv)




S_matrix_function <- function(semi_def_matrix, omega_inv){
  
  S <- matrix(data = NA, nrow = length(semi_def_matrix), ncol = length(semi_def_matrix))
  
  for (i in 1:length(semi_def_matrix)){
    for (j in 1:length(semi_def_matrix)){
      S[i,j] <- 0.5 * tr(omega_inv %*% semi_def_matrix[[i]] %*% omega_inv %*% semi_def_matrix[[j]])
    }
  }
  return(S)
}
S_matrix_function(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv)
















score_fisher_function <- function(design_matrix, semi_def_matrix, y, init_params){
  
  #-------------Initial parameters-----------------
  beta_init <- init_params[1:ncol(design_matrix)]
  sigma2_vec <- init_params[(ncol(design_matrix) + 1):length(init_params)]
  
  omega_init <- lapply(seq_along(semi_def_matrix), function(i) {            #(27.8b)
    semi_def_matrix[[i]] * sigma2_vec[i]
  })
  omega_init_inv <- chol2inv(chol(Reduce('+', omega_init)))
  
  P_init_y <- omega_init_inv %*% (y - design_matrix %*% beta_init)          #(27.17c)
  
  
  #-------------Calculating mean value parameters----------------- (27.21)
  
  M <- t(design_matrix) %*% omega_init_inv %*% design_matrix
  
  #-------------Calculating S matrix----------------- (27.22)
  
  S <- matrix(data = NA, nrow = length(semi_def_matrix), ncol = length(semi_def_matrix))
  
  for (i in 1:length(semi_def_matrix)){
    for (j in 1:length(semi_def_matrix)){
      S[i,j] <- 0.5 * tr(omega_init_inv %*% semi_def_matrix[[i]] %*% omega_init_inv %*% semi_def_matrix[[j]])
    }
  }
  
    #-------------Calculating scores-----------------
  
  score_beta <- t(design_matrix) %*% omega_init_inv %*% (y - design_matrix %*% beta_init)             # (27.10)
  
  score_sigma <- c()
  for (i in 1:length(semi_def_matrix)){
    score_sigma[i] <- -0.5 * tr(omega_init_inv %*% semi_def_matrix[[i]]) + 0.5 * t(P_init_y) %*% semi_def_matrix[[i]] %*% P_init_y
  }
  
  score <- c(score_beta, score_sigma)
  
  return(list('M' = M, 'S' = S, 'score' = score))
  
  #delta_step <- fisher_info_inv %*% score
  #return(delta_step)
}

#out <- score_fisher_function(init_params = init_params,design_matrix = design_matrix, semi_def_matrix = semi_def_matrices, y = y)

init_params <- c(1,1,1)

out <- Map(score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_params))

M_sum <- Reduce('+',lapply(out, function(x) x$M))
S_sum <- Reduce('+',lapply(out, function(x) x$S))
fisher_inv <- as.matrix(bdiag(solve(M_sum), solve(S_sum)))
score <- rowSums(sapply(out, function(x) x$score))
delta_step <- fisher_info_inv_sum %*% score
init_params <- init_params + delta_step
init_params

max_iter <- 10000000
tolerance <- 1e-6

for (iter in 1:max_iter) {
  out <- Map(score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_params))
  
  M_sum <- Reduce('+',lapply(out, function(x) x$M))
  S_sum <- Reduce('+',lapply(out, function(x) x$S))
  fisher_inv <- as.matrix(bdiag(solve(M_sum), solve(S_sum)))
  score <- rowSums(sapply(out, function(x) x$score))
  delta_step <- fisher_inv %*% score
  
  
  # Check for convergence
  if (sum((delta_step)^2) < tolerance) {
    break
  }
  
  # Update parameters for the next iteration
  init_params <- init_params + delta_step
}










tau00 <- semi_def_matrices[[1]][[1]]
tau01 <- semi_def_matrices[[1]][[2]]

tau10 <- semi_def_matrices[[2]][[1]]
tau11 <- semi_def_matrices[[2]][[2]]

beta_hat <- 1
sigma_0 <- 1
sigma_1 <- 1

V <- bdiag(sigma_0 * tau00 + sigma_1 * tau01,sigma_0 * tau10 + sigma_1 * tau11)
V0 <- bdiag(tau00, tau10)
V1 <- bdiag(tau01, tau11)
X <- matrix(rep(1,11))
y <- DF$y

score_beta <- t(X) %*% solve(V) %*% (y - X %*% beta_hat)
score_sigma_0 <- - 0.5* tr(as.matrix(solve(V) %*% V0)) + 0.5 * t(y - X %*% beta_hat) %*% solve(V) %*% V0 %*% solve(V) %*% (y - X %*% beta_hat)
score_sigma_1 <- - 0.5* tr(as.matrix(solve(V) %*% V1)) + 0.5 * t(y - X %*% beta_hat) %*% solve(V) %*% V1 %*% solve(V) %*% (y - X %*% beta_hat)

S00 <- 0.5 * tr(as.matrix(solve(V) %*% V0 %*% solve(V) %*% V0))
S10 <- 0.5 * tr(as.matrix(solve(V) %*% V0 %*% solve(V) %*% V1))
S01 <- 0.5 * tr(as.matrix(solve(V) %*% V1 %*% solve(V) %*% V0))
S11 <- 0.5 * tr(as.matrix(solve(V) %*% V1 %*% solve(V) %*% V1))

S <- matrix(c(S00, S10, S01, S11), nrow = 2,ncol = 2)

bdiag(solve(t(X) %*% solve(V) %*% X),solve(S))


tr(as.matrix(solve(V) %*% V0 %*% solve(V) %*% V0))
as.matrix(solve(V))


V_0 <- matrix(c(1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1), ncol = 4, nrow = 4, byrow = T)
V_1 <- diag(1, nrow = 4, ncol = 4)

V <- V_0 * sigma_0 + V_1 * sigma_1
V_inv <- solve (V)

S_0 <- -0.5*tr(V_inv * V_0) + 1/2 * t(y - X*beta_hat) %*% V_inv %*% V_0 %*% V_inv %*% (y - X %*% beta_hat)
S_1 <- -0.5*tr(V_inv * V_1) + 1/2 * t(y - X*beta_hat) %*% V_inv %*% V_1 %*% V_inv %*% (y - X %*% beta_hat)



V_0_0 <- matrix(c(1,1,1,1), ncol = 2, nrow = 2)
V_1_0 <- diag(1, nrow = 2, ncol = 2)
V_0 <- V_0_0 * sigma_0 + V_1_0 * sigma_1
V_0_inv <- solve(V_0)
y_0 <- y[1:2]
X_0 <- matrix(X[1:2])

S_0_0 <- -0.5*tr(V_0_inv * V_0_0) + 1/2 * t(y_0 - X_0*beta_hat) %*% V_0_inv %*% V_0_0 %*% V_0_inv %*% (y_0 - X_0 %*% beta_hat)
S_1_0 <- -0.5*tr(V_0_inv * V_1_0) + 1/2 * t(y_0 - X_0*beta_hat) %*% V_0_inv %*% V_1_0 %*% V_0_inv %*% (y_0 - X_0 %*% beta_hat)



V_0_1 <- matrix(c(1,1,1,1), ncol = 2, nrow = 2)
V_1_1 <- diag(1, nrow = 2, ncol = 2)
V_1 <- V_0_0 * sigma_0 + V_1_0 * sigma_1
V_1_inv <- solve(V_0)
y_1 <- y[3:4]
X_1 <- matrix(X[3:4])

S_0_1 <- -0.5*tr(V_1_inv * V_1_0) + 1/2 * t(y_1 - X_1*beta_hat) %*% V_1_inv %*% V_0_1 %*% V_1_inv %*% (y_1 - X_1 %*% beta_hat)
S_1_1 <- -0.5*tr(V_1_inv * V_1_1) + 1/2 * t(y_1 - X_1*beta_hat) %*% V_1_inv %*% V_1_1 %*% V_1_inv %*% (y_1 - X_1 %*% beta_hat)














S_11 <- 0.5*tr(omega_init_inv %*% gamma1 %*% omega_init_inv %*% gamma1)
#S_12 <- 0.5*tr(omega_init_inv %*% gamma1 %*% omega_init_inv %*% gamma2)
#S_21 <- 0.5*tr(omega_init_inv %*% gamma2 %*% omega_init_inv %*% gamma1)
#S_22 <- 0.5*tr(omega_init_inv %*% gamma2 %*% omega_init_inv %*% gamma2)

S <- matrix(c(S_11), nrow = 1, ncol = 1)
S_inv <- chol2inv(chol(S))


#-------------Calculating fisher information matrix----------------- (27.23)

fisher_info_inv <- as.matrix(bdiag(beta_fisher_inv, S_inv))


#-------------Calculating scores-----------------

score_beta <- t(X) %*% omega_init_inv %*% (y - X %*% beta_init)                                               # (27.10)
score_sigma_0 <- -0.5 * tr(omega_init_inv %*% gamma1) + 0.5 * t(y) %*% P_init %*% gamma1 %*% P_init %*% y     # (27.14b)
#score_sigma_1 <- -0.5 * tr(omega_init_inv %*% gamma2) + 0.5 * t(y) %*% P_init %*% gamma2 %*% P_init %*% y     # (27.14b)
score <- c(score_beta, score_sigma_0)#, score_sigma_1)


#-------------Updating parameters-----------------
init_params <- init_params + fisher_info_inv %*% score
init_params

init_params[c(2,3)] <- abs(init_params[c(2,3)])
init_params








###REML###
X <- design_matrices[[1]]
y <- outcome_list[[1]]
gamma1 <- semi_def_matrices[[1]][[1]]
gamma2 <- semi_def_matrices[[1]][[2]]

init_params <- c(1,1,1)
beta_init <- init_params[1]
omega_init <- gamma1 * init_params[2] + gamma2 * init_params[3]
omega_init_inv <- solve(omega_init)
P_init <- omega_init_inv - omega_init_inv %*% X %*% solve(t(X) %*% omega_init_inv %*% X) %*% t(X) %*% omega_init_inv


S_11 <- 0.5*tr(P_init %*% gamma1 %*% P_init %*% gamma1)
S_12 <- 0.5*tr(P_init %*% gamma1 %*% P_init %*% gamma2)
S_21 <- 0.5*tr(P_init %*% gamma2 %*% P_init %*% gamma1)
S_22 <- 0.5*tr(P_init %*% gamma2 %*% P_init %*% gamma2)

S <- matrix(c(S_11, S_12, S_21, S_22), nrow = 2, ncol = 2)
S_inv <- solve(S)

beta_fisher <- solve(t(X) %*% omega_init_inv %*% X)


fisher_info_inv <- as.matrix(bdiag(beta_fisher, S_inv))


score_beta <- t(X) %*% P_init %*% (y - X %*% beta_init)
score_sigma_0 <- -0.5 * tr(P_init %*% gamma1) + 0.5 * t(y) %*% P_init %*% gamma1 %*% P_init %*% y
score_sigma_1 <- -0.5 * tr(P_init %*% gamma2) + 0.5 * t(y) %*% P_init %*% gamma2 %*% P_init %*% y

score <- c(score_beta, score_sigma_0, score_sigma_1)


init_params <- init_params + fisher_info_inv %*% score
init_params


































#--------------Test af funktion-----------------
n_i <- nrow(matrix(model_matrix[c(1,2),]))

#Calculating covariance matrix
#1) Multiplying each semi-defnite matrix with sigma
omega_temp <- lapply(seq_along(semi_def_matrices[[1]]), function(i) {
  semi_def_matrices[[1]][[i]] * sigma2_vec[[i]]
})

#2) Summing each term in omega_temp
omega <- Reduce('+', omega_temp)


#Calculating inverse covariance matrix
omega_inv <- solve(omega)   #chol2inv(chol(omega))

#print(paste0('det(omega): ', det(omega),', log(det(omega)): ',log(det(omega))))

#Calculating log-likelihood
res <- -n_i * log(2 * pi) - 1 * log(det(omega)) - 1/2 * t(DF$y[c(1,2)] - matrix(model_matrix[c(1,2),]) %*% param_vec) %*% omega_inv %*% (DF$y[c(1,2)] - matrix(model_matrix[c(1,2),]) %*% param_vec)































data <- rnorm(8, mean = 0, sd = 1)
covariance_matrix <- diag(8)
mean_vector <- rep(0,8)
n <- length(data)
K <- dim(covariance_matrix)[1]

ll_multivariate_normal <- -n*K/2*log(2*pi) - n/2 * log(det(covariance_matrix)) - 1/2 * t(data - mean_vector) %*% covariance_matrix %*% (data - mean_vector)

ll_univariate_normal <- -n/2 * log(2*pi) - n/2 * log(covariance_matrix[1]^2) - 1/2 * sum((data-mean_vector)^2)

###########################################
#             Debugging
###########################################




design_matrices <- list(diag(1, 5, p), diag(1, 7, p), diag(1, 9, p))
semi_def_matrices <- list(list(genPositiveDefMat(dim = 5)$Sigma, genPositiveDefMat(dim = 5)$Sigma, genPositiveDefMat(dim = 5)$Sigma, genPositiveDefMat(dim = 5)$Sigma), 
                          list(genPositiveDefMat(dim = 7)$Sigma, genPositiveDefMat(dim = 7)$Sigma, genPositiveDefMat(dim = 7)$Sigma, genPositiveDefMat(dim = 7)$Sigma),
                          list(genPositiveDefMat(dim = 9)$Sigma, genPositiveDefMat(dim = 9)$Sigma, genPositiveDefMat(dim = 9)$Sigma, genPositiveDefMat(dim = 9)$Sigma))
outcome_list <- list(rnorm(5, mean = 3, sd = 1), rnorm(7, mean = 3, sd = 1), rnorm(9, mean = 3, sd = 1))
param_vec <- rnorm(p, mean = 3, sd = 1)
sigma2_vec <- rnorm(m+1,mean = 3, sd = 1)^2


log_likelihood_general(design_matrices, semi_def_matrices, outcome_list, param_vec, sigma2_vec)


#------------------------------------------------------------
n1 <- nrow(design_matrices[[1]])
n2 <- nrow(design_matrices[[2]])
#------------------------------------------------------------
omega_temp_1 <- lapply(seq_along(semi_def_matrices[[1]]), function(i) {
  semi_def_matrices[[1]][[i]] * sigma2_vec[i]
})
omega_temp_2 <- lapply(seq_along(semi_def_matrices[[2]]), function(i) {
  semi_def_matrices[[2]][[i]] * sigma2_vec[i]
})
#------------------------------------------------------------
omega1 <- Reduce('+', omega_temp_1)
omega2 <- Reduce('+', omega_temp_2)
#------------------------------------------------------------
omega1_inv <- solve(omega1)
omega2_inv <- solve(omega2)
#------------------------------------------------------------
-n1/2 * log(2 * pi) - 1/2 * log(det(omega1)) - 1/2 * t(outcome_list[[1]] - design_matrices[[1]] %*% param_vec) %*% omega1_inv %*% (outcome_list[[1]] - design_matrices[[1]] %*% param_vec)
-n_1/2 * log(2 * pi) - 1/2 * log(det(omega1)) - 1/2 * t(outcome_list[[1]] - design_matrices[[1]] %*% param_vec) %*% omega1_inv %*% (outcome_list[[1]] - design_matrices[[1]] %*% param_vec)
#------------------------------------------------------------



log_likelihood_general(design_matrices, semi_def_matrices, outcome_list, param_vec, sigma2_vec)






n_i <- nrow(design_matrices[[1]])

#Calculating covariance matrix
omega_temp <- lapply(seq_along(semi_def_matrices[[1]]), function(i) {
  semi_def_matrices[[1]][[i]] * sigma2_vec[i]
})
omega <- Reduce('+', omega_temp)
det(omega)
#Calculating inverse covariance matrix
omega_inv <- solve(omega)

#Calculating log-likelihood
-n_i/2 * log(2 * pi) - 1/2 * log(det(omega)) - 1/2 * t(outcome_list[[1]] - design_matrices[[1]] %*% param_vec) %*% omega_inv %*% (outcome_list[[1]] - design_matrices[[1]] %*% param_vec)










