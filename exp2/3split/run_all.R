rm(list = ls())
library('doParallel')
library('doRNG')
library('foreach')
library('glmnet')
library('MASS')
library("RMT4DS")
library("VGAM")
library(pracma)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
k_index <- as.numeric(args[1])
p <- as.numeric(args[2])
n_r <- as.numeric(args[3])
n_cores <- as.numeric(args[4])
n_iter_const <- as.numeric(args[5])
B <- as.numeric(args[6]) 

set.seed(1234)
registerDoParallel(cores=n_cores)
n_k <- 100
split <- 3
n <- n_r * split
gamma <- p / n_r
rho <- 0.5  
k_all <- c(
  seq(from = 0.05, to = 2, length.out = n_k / 2)[-n_k / 2], 
  seq(from = 2, to = 10, length.out = n_k / 2 + 1)
)
k <- k_all[k_index]
lambda_1 <- k_all[k_index]
lambda_2 <- k_all[k_index]
alpha_const <- runif(p); alpha_0 <- alpha_const / sqrt(sum(alpha_const^2))
beta_const <- runif(p); beta_0 <- beta_const / sqrt(sum(beta_const^2))
theta_0 = sum(alpha_0 * beta_0)
cat("Running k[", k_index, "]=", k, "\n")
cat("Running theta0:", theta_0, "\n")
mu <- c(0, 0)      
Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)

# prepare constants
mp_samples <- rgmp(n = n_iter_const, ndf = n_r, pdim = p)
c_1 <- mean(mp_samples^2 / (mp_samples + lambda_1)^2)
d_1 <- (gamma) * mean(mp_samples / (mp_samples + lambda_1)^2)
c_2 <- mean(mp_samples^2 / (mp_samples + lambda_2)^2)
d_2 <- (gamma) * mean(mp_samples / (mp_samples + lambda_2)^2)
c_3 <- mean(mp_samples / (mp_samples + lambda_1)) * mean(mp_samples / (mp_samples + lambda_2))
MP_1_b <- MP_1 <- mean(1 / (mp_samples + lambda_1))
MP_2_b <-MP_2 <- mean(1 / (mp_samples + lambda_2))
g_lambda_b <- g_lambda <- (1 - lambda_1 * MP_1_b) * (1 - lambda_2 * MP_2_b)

# Generate one dataset
X <- matrix(rnorm(n * p), n, p)
err <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
Y <- X %*% beta_0 + err[, 1]
A <- X %*% alpha_0 + err[, 2]
n1 <- floor(n / 3)
n2 <- floor(n / 3)
n3 <- n - n1 - n2
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
fit_A <- glmnet(x = X1, y = A1, family = "gaussian", alpha = 0, lambda = lambda_1 * sd_A, 
                standardize = FALSE, intercept = FALSE, thresh = 1e-25)
fit_Y <- glmnet(x = X2, y = Y2, family = "gaussian", alpha = 0, lambda = lambda_2 * sd_Y, 
                standardize = FALSE, intercept = FALSE, thresh = 1e-25)
alpha_hat <- as.vector(coef(fit_A, s = lambda_1 * sd_A))[-1]
beta_hat <- as.vector(coef(fit_Y, s = lambda_2 * sd_Y))[-1]
theta_hat <- sum(alpha_hat * beta_hat)
int_hat <- mean(A3 * Y3) - theta_hat
theta_bc_hat = theta_hat/g_lambda
pred_A <- predict(fit_A, newx = X3, s = lambda_1 * sd_A)
pred_Y <- predict(fit_Y, newx = X3, s = lambda_2 * sd_Y)
nr_hat <- mean(A3 * (Y3 - pred_Y))
if_hat <- mean((Y3 - pred_Y) * (A3 - pred_A))

int_bc_hat <- int_hat - theta_hat * (1 - g_lambda) / g_lambda
nr_bc_hat <- nr_hat - theta_hat * (lambda_2 * MP_2) / g_lambda
if_bc_hat <- if_hat - theta_hat * (lambda_1 * MP_1 * lambda_2 * MP_2) / g_lambda


generate_orthogonal_vector <- function(alpha_hat) {
  z <- rnorm(p)
  z_orth <- z - (sum(z * alpha_hat) / sum(alpha_hat^2)) * alpha_hat
  z_orth <- z_orth / sqrt(sum(z_orth^2))
  return(z_orth)
}
alpha_norm <- Norm(alpha_hat,2)
beta_norm <- Norm(beta_hat,2)
alpha_beta_dot <- sum(alpha_hat * beta_hat)
tilde_alpha <- (alpha_hat / alpha_norm) * sqrt((alpha_norm^2 - d_1) / c_1)
first_term <- (alpha_hat / alpha_norm) * (alpha_beta_dot / c_3) * (1 / sqrt((alpha_norm^2 - d_1) / c_1))
second_term_magnitude_sq <- (beta_norm^2 - d_2) / c_2 - (alpha_beta_dot / c_3 * (1 / sqrt((alpha_norm^2 - d_1) / c_1)))^2
if (second_term_magnitude_sq < 0) {
  second_term_magnitude_sq <- 0
  cat("Warning: Negative magnitude in tilde_beta computation set to 0 for k_index =", k_index, "\n")
}
z_orth <- generate_orthogonal_vector(alpha_hat)
second_term <- z_orth * sqrt(second_term_magnitude_sq)
tilde_beta <- first_term + second_term

#now we have transformed parameters tilde_alpha and tilde_beta and original estimates

run_bootstrap_iteration <- function(talpha, tbeta, rho_hat,name) {
  print(rho_hat)
  print(name)
  X_b <- matrix(rnorm(n * p), n, p)
  Sigma <- matrix(c(1, rho_hat, rho_hat, 1), nrow = 2)
  errors_b <- mvrnorm(n = n, mu = c(0, 0), Sigma = matrix(c(1, rho_hat, rho_hat, 1), nrow = 2))
  if (det(Sigma) <= 0) {
  cat("Warning: Initial Sigma matrix is not positive definite (det =", det(Sigma), ").\n")  
  }
  A_b <- X_b %*% talpha + errors_b[, 1]
  Y_b <- X_b %*% tbeta + errors_b[, 2]
  n1 <- floor(n / 3)
  n2 <- floor(n / 3)
  n3 <- n - n1 - n2
  X1_b <- X_b[1:n1, ]
  X2_b <- X_b[(n1 + 1):(n1 + n2), ]
  X3_b <- X_b[(n1 + n2 + 1):n, ]
  Y1_b <- Y_b[1:n1]
  Y2_b <- Y_b[(n1 + 1):(n1 + n2)]
  Y3_b <- Y_b[(n1 + n2 + 1):n]
  A1_b <- A_b[1:n1]
  A2_b <- A_b[(n1 + 1):(n1 + n2)]
  A3_b <- A_b[(n1 + n2 + 1):n]
  sd_A_b <- sqrt(var(A1_b) * (length(A1_b) - 1) / length(A1_b))
  sd_Y_b <- sqrt(var(Y2_b) * (length(Y2_b) - 1) / length(Y2_b))
  fit_A_b <- glmnet(x = X1_b, y = A1_b, family = "gaussian", alpha = 0, lambda = lambda_1 * sd_A_b, 
                    standardize = FALSE, intercept = FALSE, thresh = 1e-25)
  fit_Y_b <- glmnet(x = X2_b, y = Y2_b, family = "gaussian", alpha = 0, lambda = lambda_2 * sd_Y_b, 
                    standardize = FALSE, intercept = FALSE, thresh = 1e-25)
  alpha_hat_b <- as.vector(coef(fit_A_b, s = lambda_1 * sd_A_b))[-1]
  beta_hat_b <- as.vector(coef(fit_Y_b, s = lambda_2 * sd_Y_b))[-1]
  pred_A_b <- predict(fit_A_b, newx = X3_b, s = lambda_1 * sd_A_b)
  pred_Y_b <- predict(fit_Y_b, newx = X3_b, s = lambda_2 * sd_Y_b)
  theta_hat_b <- sum(alpha_hat_b * beta_hat_b)
  int_hat_b <- mean(A3_b * Y3_b) - theta_hat_b
  nr_hat_b <- mean(A3_b * (Y3_b - pred_Y_b))
  if_hat_b <- mean((Y3_b - pred_Y_b) * (A3_b - pred_A_b))
  theta_bc_hat_b <- theta_hat_b / g_lambda_b
  nr_bc_hat_b <- nr_hat_b - theta_hat_b * (lambda_2 * MP_2_b) / g_lambda_b
  if_bc_hat_b <- if_hat_b - theta_hat_b * (lambda_1 * MP_1_b * lambda_2 * MP_2_b) / g_lambda_b
  int_bc_hat_b <- int_hat_b - theta_hat_b * (1 - g_lambda_b) / g_lambda_b

  # choose which to return (BC versions as in your example)
  switch(name,
    "int" = int_bc_hat_b,
    "nr"  = nr_bc_hat_b,
    "if"  = if_bc_hat_b
  )
}



bootstrap_results <- matrix(NA, nrow = B, ncol = 3)
colnames(bootstrap_results) <- c(
  "int_bc_hat_b",
  "nr_bc_hat_b",
  "if_bc_hat_b"
)

set.seed(1234)
for (b in 1:B) {
  print(b)
  bootstrap_results[b, "int_bc_hat_b"] <- run_bootstrap_iteration(tilde_alpha, tilde_beta, int_bc_hat, name = "int")
  bootstrap_results[b, "nr_bc_hat_b"]  <- run_bootstrap_iteration(tilde_alpha, tilde_beta, nr_bc_hat, name = "nr")
  bootstrap_results[b, "if_bc_hat_b"]  <- run_bootstrap_iteration(tilde_alpha, tilde_beta, if_bc_hat, name = "if")
}


bootstrap_var_int <- var(bootstrap_results[, "int_bc_hat_b"], na.rm = TRUE)
bootstrap_var_nr <- var(bootstrap_results[, "nr_bc_hat_b"], na.rm = TRUE)
bootstrap_var_if <- var(bootstrap_results[, "if_bc_hat_b"], na.rm = TRUE)

bootstrap_df <- as.data.frame(bootstrap_results)
save(list = ls(), file = paste0("all_data_", k_index, ".RData"))
write.csv(bootstrap_df, file = paste0("bootstrap_results_k_", k_index, ".csv"), row.names = FALSE)
