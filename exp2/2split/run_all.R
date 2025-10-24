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
split <- 2
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
mu <- c(0, 0)      
Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
load("constants.RData")

gammas_k <- gammas[k_index]
zetas_k <- zetas[k_index]
nus_k <- nus[k_index]
kappas_k <- kappas[k_index]

# prepare constants
mp_samples <- rgmp(n = min(n_iter_const, 5000), ndf = n_r, pdim = p)
c1 <- mean(mp_samples^2 / (mp_samples + lambda_1)^2)
d1 <- (gamma) * mean(mp_samples / (mp_samples + lambda_1)^2)
c2 <- mean(mp_samples^2 / (mp_samples + lambda_2)^2)
d2 <- (gamma) * mean(mp_samples / (mp_samples + lambda_2)^2)
g1 <- mean(mp_samples^2 / ((mp_samples + lambda_1) * (mp_samples + lambda_2)))
g2 <- (gamma) * mean(mp_samples / ((mp_samples + lambda_1) * (mp_samples + lambda_2)))


build_transformed_coeffs <- function(alpha_hat, beta_hat, rho_hat) {
  alpha_norm <- sqrt(sum(alpha_hat^2))
  beta_norm <- sqrt(sum(beta_hat^2))
  alpha_beta_dot <- sum(alpha_hat * beta_hat)
  talpha <- (alpha_hat / alpha_norm) * sqrt((alpha_norm^2 - d1) / c1)
  proj_coeff <- ((alpha_beta_dot - g2 * rho_hat) / g1) * (1 / sqrt((alpha_norm^2 - d1) / c1))
  first_term <- (alpha_hat / alpha_norm) * proj_coeff
  second_term_magnitude_sq <- (beta_norm^2 - d2) / c2 - proj_coeff^2
  if (second_term_magnitude_sq < 0) second_term_magnitude_sq <- 0
  z <- rnorm(length(alpha_hat))
  z_orth <- z - (sum(z * alpha_hat) / sum(alpha_hat^2)) * alpha_hat
  z_orth <- z_orth / sqrt(sum(z_orth^2))
  tbeta <- first_term + z_orth * sqrt(second_term_magnitude_sq)
  return(list(talpha = talpha, tbeta = tbeta))
}

# Generate one dataset (2-split)
X <- matrix(rnorm(n * p), n, p)
err <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
Y <- X %*% beta_0 + err[, 1]
A <- X %*% alpha_0 + err[, 2]
X_tr <- X[1:n_r, ]; X_te <- X[(n_r + 1):n, ]
Y_tr <- Y[1:n_r]; Y_te <- Y[(n_r + 1):n]
A_tr <- A[1:n_r]; A_te <- A[(n_r + 1):n]
sd_A <- sqrt(var(A_tr) * (n_r - 1) / n_r)
sd_Y <- sqrt(var(Y_tr) * (n_r - 1) / n_r)
fit_A <- glmnet(x = X_tr, y = A_tr, family = "gaussian", alpha = 0, lambda = lambda_1 * sd_A, 
                standardize = FALSE, intercept = FALSE, thresh = 1e-25)
fit_Y <- glmnet(x = X_tr, y = Y_tr, family = "gaussian", alpha = 0, lambda = lambda_2 * sd_Y, 
                standardize = FALSE, intercept = FALSE, thresh = 1e-25)
alpha_hat <- as.vector(coef(fit_A, s = lambda_1 * sd_A))[-1]
beta_hat <- as.vector(coef(fit_Y, s = lambda_2 * sd_Y))[-1]
theta_hat <- sum(alpha_hat * beta_hat)
pred_A <- predict(fit_A, newx = X_te, s = lambda_1 * sd_A)
pred_Y <- predict(fit_Y, newx = X_te, s = lambda_2 * sd_Y)
 int_hat <- mean(A_te * Y_te) - theta_hat
nr_hat <- mean(A_te * (Y_te - pred_Y))
if_hat <- mean((Y_te - pred_Y) * (A_te - pred_A))
rho_bc_hat_NR <- (nr_hat - theta_hat * (1 - nus_k) / gammas_k) * (1 / (1 - zetas_k * (1 - nus_k) / gammas_k))
rho_bc_hat_DR <- (if_hat - theta_hat * (kappas_k) / gammas_k) * (1 / (1 + zetas_k - zetas_k * kappas_k / gammas_k))

 # Transformed coefficients per estimator
 tcoef_INT <- build_transformed_coeffs(alpha_hat, beta_hat, rho_bc_hat_INT)
tcoef_NR  <- build_transformed_coeffs(alpha_hat, beta_hat, rho_bc_hat_NR)
tcoef_DR  <- build_transformed_coeffs(alpha_hat, beta_hat, rho_bc_hat_DR)


run_bootstrap_iteration <- function(talpha, tbeta, rho_hat,name) {
  print(rho_hat)
  print(name)
  Sigma <- matrix(c(1, rho_hat, rho_hat, 1), nrow = 2)

  if (det(Sigma) <= 0) {
  cat("Warning: Initial Sigma matrix is not positive definite (det =", det(Sigma), ").\n")  
  }
    Xb <- matrix(rnorm(n * p), n, p)
    Eb <- mvrnorm(n = n, mu = c(0, 0), Sigma = matrix(c(1, rho_hat, rho_hat, 1), nrow = 2))
    Ab <- Xb %*% talpha + Eb[, 2]
    Yb <- Xb %*% tbeta + Eb[, 1]
    Xtr <- Xb[1:n_r, ]; Xte <- Xb[(n_r + 1):n, ]
    Atr <- Ab[1:n_r]; Ate <- Ab[(n_r + 1):n]
    Ytr <- Yb[1:n_r]; Yte <- Yb[(n_r + 1):n]
    sd_Atr <- sqrt(var(Atr) * (n_r - 1) / n_r)
    sd_Ytr <- sqrt(var(Ytr) * (n_r - 1) / n_r)
    lamA <- k * sd_Atr
    lamY <- k * sd_Ytr
    fitA <- glmnet(x = Xtr, y = Atr, family = "gaussian", alpha = 0, lambda = lamA,
                   standardize = FALSE, intercept = FALSE, thresh = 1e-25)
    fitY <- glmnet(x = Xtr, y = Ytr, family = "gaussian", alpha = 0, lambda = lamY,
                   standardize = FALSE, intercept = FALSE, thresh = 1e-25)
    alpha_hat_b <- as.vector(coef(fitA, s = lamA))[-1]
    beta_hat_b <- as.vector(coef(fitY, s = lamY))[-1]
    pred_Ate <- predict(fitA, newx = Xte, s = lamA)
    pred_Yte <- predict(fitY, newx = Xte, s = lamY)
    theta_hat_b <- sum(alpha_hat_b * beta_hat_b)
    int_b <- mean(Ate * Yte) - theta_hat_b
    nr_b <- mean(Ate * (Yte - pred_Yte))
    dr_b <- mean((Yte - pred_Yte) * (Ate - pred_Ate))

  switch(name,
  "int" = (int_b - theta_hat_b * (1 - gammas_k) / gammas_k) * (gammas_k / (gammas_k - zetas_k)),
  "nr"  = (nr_b - theta_hat_b * (1 - nus_k) / gammas_k) * (1 / (1 - zetas_k * (1 - nus_k) / gammas_k)),
  "if"  = (dr_b - theta_hat_b * (kappas_k) / gammas_k) * (1 / (1 + zetas_k - zetas_k * kappas_k / gammas_k))
  )
}
bootstrap_results <- matrix(NA, nrow = B, ncol = 3)
bootstrap_results <- matrix(NA, nrow = B, ncol = 2)
colnames(bootstrap_results) <- c(
  "int_bc_hat_b",
  "nr_bc_hat_b",
  "if_bc_hat_b"
)

set.seed(1234)
for (b in 1:B) {
  print(b)
  bootstrap_results[b, "int_bc_hat_b"] <- run_bootstrap_iteration(tcoef_INT$talpha, tcoef_INT$tbeta,rho_bc_hat_INT, name = "int")
  bootstrap_results[b, "nr_bc_hat_b"]  <- run_bootstrap_iteration(tcoef_NR$talpha, tcoef_NR$tbeta, rho_bc_hat_NR, name = "nr")
  bootstrap_results[b, "if_bc_hat_b"]  <- run_bootstrap_iteration(tcoef_DR$talpha, tcoef_DR$tbeta, rho_bc_hat_DR, name = "if")
}

bootstrap_var_int <- var(bootstrap_results[, "int_bc_hat_b"], na.rm = TRUE)
bootstrap_var_nr <- var(bootstrap_results[, "nr_bc_hat_b"], na.rm = TRUE)
bootstrap_var_if <- var(bootstrap_results[, "if_bc_hat_b"], na.rm = TRUE)

bootstrap_df <- as.data.frame(bootstrap_results)
save(list = ls(), file = paste0("all_data_", k_index, ".RData"))
write.csv(bootstrap_df, file = paste0("bootstrap_results_k_", k_index, ".csv"), row.names = FALSE)
