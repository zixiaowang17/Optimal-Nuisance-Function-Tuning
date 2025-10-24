rm(list = ls())

library('doParallel')
library('doRNG')
library('foreach')
library('glmnet')
library('MASS')

source('settings.R')

registerDoParallel(cores=n_cores)

alpha <- alpha_const; beta <- beta_const
n_iter <- n_iter_const


################################################################################
## Get taus
################################################################################
set.seed(1234)
cat("Simulating taus...\n")
sim_res_taus <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  e_vals <- eigen(t(X) %*% X / n_r, symmetric = TRUE, only.values = TRUE)$values
  
  for (j in 1:n_lambda){
    res[1, j] <- mean(e_vals^2  / (e_vals + k_all[j])^4)
  }
  return(data.frame(res))
}
taus <- unname(colMeans(sim_res_taus) * p / n_r)


################################################################################
## Get iotas
################################################################################
set.seed(1234)
cat("Simulating iotas...\n")
sim_res_iotas <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Sigma <- t(X) %*% X / n_r
  
  LHS <- t(alpha) %*% Sigma
  for (j in 1:n_lambda){
    Center <- Sigma + (k_all[j] * diag(p))
    Center_temp <- solve(Center %*% Center, Sigma)
    res[1, j] <- LHS %*% Center_temp %*% Center_temp %*% alpha
  }
  return(data.frame(res))
}
iotas <- unname(colMeans(sim_res_iotas) / sum(alpha * alpha))



################################################################################
## Get Ts
################################################################################
cat("Simulating Ts...\n")
set.seed(1234)
sim_res_T1s <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Sigma <- t(X) %*% X / n_r
  
  LHS <- t(alpha) %*% Sigma
  RHS <- Sigma %*% beta
  for (j in 1:n_lambda){
    Center <- Sigma + (k_all[j] * diag(p))
    res[1, j] <- LHS %*% solve(Center %*% Center, RHS)
  }
  return(data.frame(res))
}
T1s <- apply(sim_res_T1s, 2, var)


set.seed(1234)
cat("Simulating T4s...\n")
sim_res_T4s <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  eps <- mvrnorm(n = n_r, mu = c(0, 0), 
                 Sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  
  Sigma <- t(X) %*% X / n_r
  
  LHS <- t(t(X) %*% eps[, 1]) / n_r
  RHS <- t(X) %*% eps[, 2] / n_r
  for (j in 1:n_lambda){
    Center <- Sigma + (k_all[j] * diag(p))
    res[1, j] <- LHS %*% solve(Center %*% Center, RHS)
  }
  return(data.frame(res))
}
T4s <- apply(sim_res_T4s, 2, var)


save(taus, iotas, T1s, T4s, file = 'constants2.RData')
