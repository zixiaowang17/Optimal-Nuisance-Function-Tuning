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
## Get zetas
################################################################################
cat("Simulating zetas...\n")
set.seed(1234)
sim_res_zetas <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  e_vals <- eigen(cov(X), symmetric = TRUE, only.values = TRUE)$values
  
  for (j in 1:n_lambda){
    res[1, j] <- mean(e_vals  / (e_vals + k_all[j])^2)
  }
  return(data.frame(res))
}
zetas <- unname(colMeans(sim_res_zetas) * p / n_r)


################################################################################
## Get gammas
################################################################################
cat("Simulating gammas...\n")
set.seed(1234)
sim_res_gammas <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Sigma <- cov(X)
  
  LHS <- t(alpha) %*% Sigma
  RHS <- Sigma %*% beta
  for (j in 1:n_lambda){
    Center <- Sigma + (k_all[j] * diag(p))
    res[1, j] <- LHS %*% solve(Center %*% Center, RHS)
  }
  return(data.frame(res))
}
gammas <- unname(colMeans(sim_res_gammas) / sum(alpha * beta))


################################################################################
## Get nus
################################################################################
cat("Simulating nus...\n")
sim_res_nus <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Y <- X %*% beta + rnorm(n_r)
  
  sd_Y <- sqrt(var(Y) * (n_r - 1) / n_r)[1, 1]
  lambda_all <- k_all * n_r
  glmnet_lambda_Y <- lambda_all * sd_Y / n_r
  
  fit_Y <- glmnet(x = X, y = Y, family = "gaussian", alpha = 0, lambda = glmnet_lambda_Y, 
                  standardize = FALSE, intercept = F, thresh = 1e-25)
  
  X_new <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Sigma <- cov(X_new)
  
  for (j in 1:n_lambda){
    beta_Y <- as.vector(coef(fit_Y, s = glmnet_lambda_Y[j]))[-1]
    res[1, j] <- sum(alpha * (Sigma %*% beta_Y))
  }
  return(data.frame(res))
}
nus <- unname(colMeans(sim_res_nus) / sum(alpha * beta))


################################################################################
## Get kappas
################################################################################
cat("Simulating kappas...\n")
sim_res_kappas <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Sigma <- cov(X)
  
  for (j in 1:n_lambda){
    k_I <- diag(k_all[j], nrow = p, ncol = p)
    RHS <- Sigma %*% solve(Sigma + k_I, alpha) - alpha
    temp <- Sigma %*% solve(Sigma + k_I, RHS) - RHS
    res[1, j] <- sum(beta * temp)
  }
  return(data.frame(res))
}
kappas <- unname(colMeans(sim_res_kappas) / sum(alpha * beta))


################################################################################
## Get gammas prime
################################################################################
cat("Simulating gammas prime...\n")
sim_res_gammas_p <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  res <- matrix(NA, nrow = 1, ncol = n_lambda)
  
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Sigma <- cov(X)
  
  RHS <- Sigma %*% alpha
  for (j in 1:n_lambda){
    res[1, j] <- t(alpha) %*% solve(Sigma + (k_all[j] * diag(p)), RHS)
  }
  return(data.frame(res))
}

gammas_p <- unname(colMeans(sim_res_gammas_p) / sum(alpha * alpha))


save(zetas, gammas, nus, kappas, gammas_p, file = 'constants.RData')
