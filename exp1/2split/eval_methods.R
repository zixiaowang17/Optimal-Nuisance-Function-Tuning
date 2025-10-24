rm(list = ls())

library('doParallel')
library('doRNG')
library('foreach')
library('glmnet')
library('MASS')
source('settings.R')

load('constants.RData')

set.seed(1234)
registerDoParallel(cores=n_cores)

alpha <- alpha_sim; beta <- beta_sim
n_iter <- n_iter_sim

sim_res <- foreach(i = 1:n_iter) %dorng% {
  res_int <- res_int_db <- res_nr <- res_nr_db <- res_if <- res_if_db <- matrix(NA, nrow = 1, ncol = n_lambda)
  res_ols <- NA
  
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
  
  X_test <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
  eps_test <- mvrnorm(n = n_r, mu = c(0, 0), Sigma = Sigma)
  Y_test <- X_test %*% beta + eps_test[, 1]
  A_test <- X_test %*% alpha  + eps_test[, 2]
  
  for (j in 1:n_lambda){
    beta_Y <- as.vector(coef(fit_Y, s = glmnet_lambda_Y[j]))[-1]
    beta_A <- as.vector(coef(fit_A, s = glmnet_lambda_A[j]))[-1]
    pred_Y <- predict(fit_Y, newx = X_test, s = glmnet_lambda_Y[j])
    pred_A <- predict(fit_A, newx = X_test, s = glmnet_lambda_A[j])
    
    # Integral-based plug-in estimator
    int <- mean(A_test * Y_test) - sum(beta_Y * beta_A)
    mult_int <- (gammas[j] / (gammas[j] - zetas[j]))
    res_int[1, j] <- int
    res_int_db[1, j] <- (int - sum(beta_Y * beta_A) * (1 - gammas[j]) / gammas[j]) * mult_int
    
    # Newey-Robins plug-in estimator
    est_nr <- mean(A_test * (Y_test - pred_Y))
    mult_nr <- 1 / (1 - zetas[j] * (1 - nus[j]) / gammas[j])
    res_nr[1, j] <- est_nr
    res_nr_db[1, j] <- (est_nr - sum(beta_Y * beta_A) * (1 - nus[j]) / gammas[j]) * mult_nr
    
    # First-order estimator
    est_if <- mean((Y_test - pred_Y) * (A_test - pred_A))
    mult_if <- 1 / (1 + zetas[j] - zetas[j] * kappas[j] / gammas[j])
    res_if[1, j] <- est_if
    res_if_db[1, j] <- (est_if - sum(beta_Y * beta_A) * kappas[j] / gammas[j]) * mult_if
  }
  
  if (include_ols){
    X_test2 <- matrix(rnorm(n_r * p), nrow = n_r, ncol = p)
    Y_test2 <- X_test2 %*% beta + rnorm(n_r)
    
    fit_A_ols <- lm(A ~ X - 1)
    fit_Y_ols <- lm(Y_test2 ~ X_test2 - 1)
    res_ols <- mean(A_test * Y_test) - sum(coef(fit_Y_ols) * coef(fit_A_ols))
  }
  return(list(res_int = res_int, res_int_db = res_int_db, 
              res_if = res_if, res_if_db = res_if_db, 
              res_nr = res_nr, res_nr_db = res_nr_db, 
              res_ols = res_ols))
}

# Consolidate results
res_int <- res_int_db <- res_nr <- res_nr_db <- res_if <- res_if_db <- matrix(NA, nrow = n_iter, ncol = n_lambda)
res_ols <- matrix(NA, nrow = n_iter, ncol = 1)
for (i in 1:n_iter){
  res_int[i, ] <- sim_res[[i]]$res_int
  res_int_db[i, ] <- sim_res[[i]]$res_int_db
  
  res_nr[i, ] <- sim_res[[i]]$res_nr
  res_nr_db[i, ] <- sim_res[[i]]$res_nr_db
  
  res_if[i, ] <- sim_res[[i]]$res_if
  res_if_db[i, ] <- sim_res[[i]]$res_if_db
  
  res_ols[i, ] <- sim_res[[i]]$res_ols
}

save(res_int, res_int_db, res_nr, res_nr_db, res_if, res_if_db, res_ols,
     file = 'Sim-Res.RData')
