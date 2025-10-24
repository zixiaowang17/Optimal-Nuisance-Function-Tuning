args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
p <- args[2]
rho = 0.5
n_k <- 100
k_all <- c(
  seq(from = 0.05, to = 2, length.out = n_k / 2)[-n_k / 2], 
  seq(from = 2, to = 10, length.out = n_k / 2 + 1)
)
lambda_1 <- lambda_2 <- k_all 

# Preallocate
theta_k     <- numeric(n_k)
int_k       <- numeric(n_k)
nr_k        <- numeric(n_k)
if_k        <- numeric(n_k)
theta_bc_k  <- numeric(n_k)
int_bc_k    <- numeric(n_k)
nr_bc_k     <- numeric(n_k)
if_bc_k     <- numeric(n_k)
theta_0     <- numeric(n_k)
theta_k_var    <- numeric(n_k)
int_k_var      <- numeric(n_k)
nr_k_var       <- numeric(n_k)
if_k_var       <- numeric(n_k)
theta_bc_k_var <- numeric(n_k)
int_bc_k_var   <- numeric(n_k)
nr_bc_k_var    <- numeric(n_k)
if_bc_k_var    <- numeric(n_k)
MP_1_k <- numeric(n_k)
MP_2_k <- numeric(n_k)
g_lambda_k <- numeric(n_k)

# ---- Read files and compute summary ----
for (i in 1:n_k) {  
  file_path <- file.path(dir_path, paste0("result_k_", i, ".csv"))

  if (file.exists(file_path)) {
    df <- read.csv(file_path)
    theta_k[i]     <- mean(df$theta_hat)
    int_k[i]       <- mean(df$int_hat)
    nr_k[i]        <- mean(df$nr_hat)
    if_k[i]        <- mean(df$if_hat)
    theta_bc_k[i]  <- mean(df$theta_bc_hat)
    int_bc_k[i]    <- mean(df$int_bc_hat)
    nr_bc_k[i]     <- mean(df$nr_bc_hat)
    if_bc_k[i]     <- mean(df$if_bc_hat)
    theta_0[i] <- mean(df$theta_0)
    MP_1_k[i] <- mean(df$MP_1)
    MP_2_k[i] <- mean(df$MP_2)
    g_lambda_k[i] <- mean(df$g_lambda)
    theta_k_var[i]    <- var(df$theta_hat)
    int_k_var[i]      <- var(df$int_hat)
    nr_k_var[i]       <- var(df$nr_hat)
    if_k_var[i]       <- var(df$if_hat)
    theta_bc_k_var[i] <- var(df$theta_bc_hat)
    int_bc_k_var[i]   <- var(df$int_bc_hat)
    nr_bc_k_var[i]    <- var(df$nr_bc_hat)
    if_bc_k_var[i]    <- var(df$if_bc_hat)
  } else {
    warning(paste("Missing file for k =", i))
  }
}



pdf(paste0("3_p_",p,'_Naive-vs-Debiased.pdf'), height = 4, width = 12)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2),    # inner margins
    oma = c(0.5, 0.5, 0.5, 0.5),  # outer margins
    mgp = c(3.4, 1.1, 0))
cex_main <- 1.5
cex_lab  <- 1.5
cex_axis <- 1.35

plot(k_all, int_k - rho, ylim = c(-1, 1), ylab = 'Bias', xlab = expression(lambda),
      pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
      main = 'Integral-based Estimator (Three Splits)')
 points(k_all, int_bc_k - rho, col = 'blue', pch = 19, cex = 0.75)
 lines(k_all, theta_0 * (1-g_lambda_k) , col = 'red')
 lines(k_all, rep(0, n_k), col = 'red')

 plot(k_all, nr_k - rho, ylim = c(-1, 1), ylab = 'Bias', xlab = expression(lambda),
      pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
      main = 'Newey and Robins Estimator (Three Splits)')
 points(k_all, nr_bc_k - rho, col = 'blue', pch = 19, cex = 0.75)
 lines(k_all, theta_0 * (lambda_2*MP_2_k), col = 'red')
 lines(k_all, rep(0, n_k), col = 'red')

plot(k_all, if_k - rho, ylim = c(-1, 1), ylab = 'Bias', xlab = expression(lambda),
     pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     main = 'Doubly Robust Estimator (Three Splits)')
points(k_all, if_bc_k - rho, col = 'blue', pch = 19, cex = 0.75)
lines(k_all, theta_0  *lambda_1*MP_1_k*lambda_2*MP_2_k, col = 'red')
lines(k_all, rep(0, n_k), col = 'red')


load('Sim-Res-Pred.RData')

mean_pred_p <- colMeans(res_p)
mean_pred_b <- colMeans(res_b)


pdf(paste0("3_p_",p,"_Variances.pdf"), height = 4, width = 12)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2),    # inner margins
    oma = c(0.5, 0.5, 0.5, 0.5),  # outer margins
    mgp = c(3.4, 1.1, 0))
cex_main <- 1.5
cex_lab  <- 1.5
cex_axis <- 1.35

plot(k_all, int_bc_k_var, ylab = 'Variance', xlab = expression(lambda),
     main =  'Integral-based Estimator (Three Splits)', cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75)
abline(v = k_all[which.min(int_bc_k_var)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

plot(k_all, nr_bc_k_var, ylab = 'Variance', xlab = expression(lambda),
     main = 'Newey and Robins Estimator (Three Splits)', cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75)
abline(v = k_all[which.min(nr_bc_k_var)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

plot(k_all, if_bc_k_var, ylab = 'Variance', xlab = expression(lambda),
     main = 'Doubly Robust Estimator (Three Splits)', cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75)
abline(v = k_all[which.min(if_bc_k_var)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

dev.off()

if (p == 1000) {
pdf(paste0("cut3p.pdf"), height = 4, width = 12)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2),    # inner margins
    oma = c(0.5, 0.5, 0.5, 0.5),  # outer margins
    mgp = c(3.4, 1.1, 0))
cex_main <- 1.5
cex_lab  <- 1.5
cex_axis <- 1.35

plot(k_all, int_bc_k_var, ylab = 'Variance', xlab = expression(lambda),
     main =  'Integral-based Estimator (Three Splits)',cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75,cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     ylim = c(0, 0.05), 
     xlim = c(0, 4))
abline(v = k_all[which.min(int_bc_k_var)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

plot(k_all, nr_bc_k_var, ylab = 'Variance', xlab = expression(lambda),
     main = 'Newey and Robins Estimator (Three Splits)', cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75,cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     ylim = c(0, 0.05), 
     xlim = c(0, 4))
abline(v = k_all[which.min(nr_bc_k_var)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

plot(k_all, if_bc_k_var, ylab = 'Variance', xlab = expression(lambda),
     main = 'Doubly Robust Estimator (Three Splits)', cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75,cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     ylim = c(0, 0.05), 
     xlim = c(0, 4))
abline(v = k_all[which.min(if_bc_k_var)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

dev.off()
}
