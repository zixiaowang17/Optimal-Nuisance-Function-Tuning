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
bootstrap_var_int_k     <- numeric(n_k)
bootstrap_var_nr_k     <- numeric(n_k)
bootstrap_var_if_k     <- numeric(n_k)
emperial_var_int_k     <- numeric(n_k)
emperial_var_nr_k     <- numeric(n_k)
emperial_var_if_k     <- numeric(n_k)

for (i in 1:n_k) {  
  file_path <- file.path(dir_path, paste0("bootstrap_results_k_", i, ".csv"))
  print(i)

  if (file.exists(file_path)) {
    df <- read.csv(file_path)
    bootstrap_var_int_k[i]     <- var(df$int_bc_hat_b)
    bootstrap_var_nr_k[i]     <- var(df$nr_bc_hat_b)
    bootstrap_var_if_k[i]     <- var(df$if_bc_hat_b)
  } else {
    warning(paste("Missing file for k =", i))
  }
}
for (i in 1:n_k) {  
  file_path <- file.path(paste0("result_k_", i, ".csv"))
  print(i)
  if (file.exists(file_path)) {
    df <- read.csv(file_path)
    emperial_var_int_k[i]     <- var(df$int_bc_hat)
    emperial_var_nr_k[i]     <- var(df$nr_bc_hat)
    emperial_var_if_k[i]     <- var(df$if_bc_hat)
  } else {
    warning(paste("Missing file for k =", i))
  }
}

#std of the ratio
summary(sqrt(emperial_var_int_k/bootstrap_var_int_k))
summary(sqrt(emperial_var_nr_k/bootstrap_var_nr_k))
summary(sqrt(emperial_var_if_k/bootstrap_var_if_k))
pdf(paste0("Boot_Variances_p_",p,".pdf"), height = 4, width = 12)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2),    # inner margins
    oma = c(0.5, 0.5, 0.5, 0.5),  # outer margins
    mgp = c(3.4, 1.1, 0))
cex_main <- 1.5
cex_lab  <- 1.5
cex_axis <- 1.35

p_num <- as.numeric(p)
ymax <-   max(c(
    bootstrap_var_int_k, emperial_var_int_k,
    bootstrap_var_nr_k,  emperial_var_nr_k,
    bootstrap_var_if_k,  emperial_var_if_k
  ), na.rm = TRUE)


# INT (integral-based plug-in)
plot(k_all, bootstrap_var_int_k, ylab = 'Variance', xlab = expression(lambda),
     main = 'Integral-Based Estimator (Two Splits)',
     pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     ylim = range(c(0, ymax), na.rm = TRUE))
lines(k_all, emperial_var_int_k, col = 'red', lty = 2)
# NR (Newey and Robins)
plot(k_all, bootstrap_var_nr_k, ylab = 'Variance', xlab = expression(lambda),
     main = 'Newey and Robins Estimator (Two Splits)',
     pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     ylim = range(c(0, ymax), na.rm = TRUE))
lines(k_all, emperial_var_nr_k, col = 'red', lty = 2)
# DR (doubly robust)
plot(k_all, bootstrap_var_if_k, ylab = 'Variance', xlab = expression(lambda),
     main = 'Doubly Robust Estimator (Two Splits)',
     pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     ylim = range(c(0, ymax), na.rm = TRUE))
lines(k_all, emperial_var_if_k, col = 'red', lty = 2)
dev.off()
