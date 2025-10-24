rm(list = ls())
library('doParallel')
library('doRNG')
library('foreach')
library('glmnet')
library('MASS')
library("RMT4DS")

# Parse command-line argument
args <- commandArgs(trailingOnly = TRUE)
k_index <- as.numeric(args[1])
p  <- as.numeric(args[2])
n_r <- as.numeric(args[3])
n_cores <-as.numeric(args[4])
n_iter_const <- as.numeric(args[5])
n_iter_sim <- as.numeric(args[6])

set.seed(1234)
registerDoParallel(cores=n_cores)
n_k <- 100
split = 3
n = n_r*split
gamma = p/n_r
rho <- 0.5  
 #number of different values of lambda we are fitting models for 
k_all <- c(
  seq(from = 0.05, to = 2, length.out = n_k / 2)[-n_k / 2], 
  seq(from = 2, to = 10, length.out = n_k / 2 + 1)
)
k <- k_all[k_index]
lambda_1 = k_all[k_index]
lambda_2 = k_all[k_index]
 
alpha_const <- runif(p); alpha_0 <- alpha_const / sqrt(sum(alpha_const^2))
beta_const <- runif(p); beta_0 <- beta_const / sqrt(sum(beta_const^2))
theta_0 = sum(alpha_0 * beta_0)

  cat("Running k[", k_index, "]=",k, "\n")
  cat("Running alpha0:", alpha_0,"\n")
  cat("Running beta0:", beta_0,"\n")
  cat("Running theta0:", theta_0,"\n")

a <- (1 - sqrt(gamma))^2
b <- (1 + sqrt(gamma))^2

norm = function(x)
  return(sqrt(sum(x**2)))
mu = c(0, 0)      
Sigma = matrix(c(1, rho,
                 rho, 1), nrow = 2)
theta_hat <- int_hat <- nr_hat <- if_hat <- theta_bc_hat <- int_bc_hat <- nr_bc_hat <-  if_bc_hat <- numeric(n_iter_sim)  
MP_1 <- MP_2 <- g_lambda <- numeric(n_iter_sim)  


results <- foreach(j = 1:n_iter_sim, .combine = rbind, .packages = c("MASS")) %dorng% {
  X = matrix(rnorm(n * p), n, p)
  err = mvrnorm(n = n, mu = mu, Sigma = Sigma)
  Y = X %*% beta_0 + err[,1]
  A = X %*% alpha_0 + err[,2]
  n1 <- floor(n / 3)
  n2 <- floor(n / 3)
  n3 <- n - n1 - n2  # Remaining samples
  X1 <- X[1:n1, ]
  X2 <- X[(n1 + 1):(n1 + n2), ]
  X3 <- X[(n1 + n2 + 1):n, ]
  Y1 <- Y[1:n1]
  Y2 <- Y[(n1 + 1):(n1 + n2)]
  Y3 <- Y[(n1 + n2 + 1):n]
  A1 <- A[1:n1]
  A2 <- A[(n1 + 1):(n1 + n2)]
  A3 <- A[(n1 + n2 + 1):n]  
  

  sd_A <- sqrt(var(A1) * (n_r - 1) / n_r)
  sd_Y <- sqrt(var(Y2) * (n_r - 1) / n_r)
  fit_A <- glmnet(x = X1, y = A1, family = "gaussian", alpha = 0, lambda = lambda_1 * sd_A , 
                  standardize = FALSE, intercept = F, thresh = 1e-25)
  fit_Y <- glmnet(x = X2, y = Y2, family = "gaussian", alpha = 0, lambda = lambda_2 * sd_Y , 
                  standardize = FALSE, intercept = F, thresh = 1e-25)
  pred_A <- predict(fit_A, newx = X3, s = lambda_1 * sd_A)
  pred_Y <- predict(fit_Y, newx = X3, s = lambda_2 * sd_Y)
  alpha_hat <- as.vector(coef(fit_A, s = lambda_1 * sd_A ))[-1]
  beta_hat <- as.vector(coef(fit_Y, s = lambda_2 * sd_Y))[-1]

  mp_samples <- rgmp(n = n_iter_const, ndf = n_r, pdim = p)
  MP_1[j] = mean((1/(mp_samples + lambda_1)))
  MP_2[j] = mean((1/(mp_samples + lambda_2)))
  g_lambda[j] = (1-lambda_1*MP_1[j])*(1-lambda_2*MP_2[j])

  theta_hat[j] = sum(alpha_hat * beta_hat)
  int_hat[j] = mean(A3 * Y3) - sum(alpha_hat * beta_hat)
  nr_hat[j] = mean(A3 * (Y3 - pred_Y))
  if_hat [j] = mean((Y3 - pred_Y) * (A3 - pred_A))
  theta_bc_hat[j] = theta_hat[j]/g_lambda[j]
  nr_bc_hat[j] = nr_hat[j] - theta_hat[j]*(lambda_2*MP_2[j])/g_lambda[j]
  if_bc_hat[j] = if_hat[j] -theta_hat[j]*(lambda_1*MP_1[j]*lambda_2*MP_2[j])/g_lambda[j]
  int_bc_hat[j] = int_hat[j] - theta_hat[j]*(1-g_lambda[j])/g_lambda[j]
 c(
  theta_hat[j],
  int_hat[j],
  nr_hat[j],
  if_hat[j],
  theta_bc_hat[j],
  int_bc_hat[j],
  nr_bc_hat[j],
  if_bc_hat[j],
  MP_1[j],
  MP_2[j],
  g_lambda[j]
)
}


results_df <- as.data.frame(results)
colnames(results_df) <- c("theta_hat", "int_hat", "nr_hat", "if_hat", "theta_bc_hat", "int_bc_hat", "nr_bc_hat", "if_bc_hat","MP_1","MP_2","g_lambda")
results_df$theta_0 <- theta_0
write.csv(results_df, file = paste0("result_k_", k_index, ".csv"), row.names = FALSE)