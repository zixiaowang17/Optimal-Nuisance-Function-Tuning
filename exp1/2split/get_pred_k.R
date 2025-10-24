rm(list = ls())

library('doParallel')
library('doRNG')
library('foreach')
library('glmnet')
library('MASS')

source('settings.R')

set.seed(1234)
registerDoParallel(cores=n_cores)

alpha <- alpha_sim; beta <- beta_sim
n_iter <- n_iter_sim

sim_res <- foreach(i = 1:n_iter) %dorng% {
  res_p <- res_b <- matrix(NA, nrow = 1, ncol = n_lambda)

  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  eps <- mvrnorm(n = n_r, mu = c(0, 0), Sigma = Sigma)
  Y <- X %*% beta + eps[, 1]
  A <- X %*% alpha  + eps[, 2]
  
  sd_A <- sqrt(var(A) * (n_r - 1) / n_r)[1, 1]
  sd_Y <- sqrt(var(Y) * (n_r - 1) / n_r)[1, 1]
  
  lambda_all <- k_all * n_r
  glmnet_lambda_A <- lambda_all * sd_A / n_r
  glmnet_lambda_Y <- lambda_all * sd_Y / n_r
  
  fit_Y <- glmnet(x = X, y = Y, family = "gaussian", alpha = 0, lambda = glmnet_lambda_Y, 
                  standardize = FALSE, intercept = F, thresh = 1e-25)
  fit_A <- glmnet(x = X, y = A, family = "gaussian", alpha = 0, lambda = glmnet_lambda_A, 
                  standardize = FALSE, intercept = F, thresh = 1e-25)
  
  
  # Evaluate methods
  X <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  Y_mu <- X %*% beta
  A_mu <- X %*% alpha
  
  # Estimating p for each candidate model
  for (j in 1:n_lambda){
    A_hat <- predict(fit_A, newx = X, s = glmnet_lambda_Y[j])
    res_p[1, j] <- mean((A_hat - A_mu)^2)
  }
  
  # Estimating b for each candidate model
  for (j in 1:n_lambda){
    Y_hat <- predict(fit_Y, newx = X, s = glmnet_lambda_Y[j])
    res_b[1, j] <- mean((Y_hat - Y_mu)^2)
  }

  return(list(res_p = res_p, res_b = res_b))
}


# Consolidate results
res_p <- res_b <- matrix(NA, nrow = n_iter, ncol = n_lambda)
for (i in 1:n_iter){
  res_p[i, ] <- sim_res[[i]]$res_p
  res_b[i, ] <- sim_res[[i]]$res_b
}

save(res_p, res_b, file = 'Sim-Res-Pred.RData')
