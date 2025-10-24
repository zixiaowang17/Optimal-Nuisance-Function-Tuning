n_iter_const <- 10000
n_iter_sim <- 10000
n <- 1000 
n_r <- round(n / 2)
p <- 250
rho <- 0.5
n_lambda <- 100
k_all <- c(
  seq(from = 0.05, to = 2, length.out = n_lambda / 2)[-n_lambda / 2], 
  seq(from = 2, to = 10, length.out = n_lambda / 2 + 1))
include_ols <- TRUE

set.seed(1234)
alpha_const <- runif(p); alpha_const <- alpha_const / sqrt(sum(alpha_const^2))
beta_const <- runif(p); beta_const <- beta_const / sqrt(sum(beta_const^2))

alpha_sim <- runif(p); alpha_sim <- alpha_sim / sqrt(sum(alpha_sim^2))
beta_sim <- runif(p); beta_sim <- beta_sim / sqrt(sum(beta_sim^2))

n_cores <- 110