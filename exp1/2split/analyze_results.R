rm(list = ls())
source('settings.R')
load('constants.RData')
load('constants2.RData')
load('Sim-Res.RData')
load('Sim-Res-Pred.RData')

alpha <- alpha_sim
beta <- beta_sim 

mean_int <- colMeans(res_int)
mean_int_db <- colMeans(res_int_db)
var_int <- apply(res_int, 2, var)
var_int_db <- apply(res_int_db, 2, var)

mean_nr <- colMeans(res_nr)
mean_nr_db <- colMeans(res_nr_db)
var_nr_db <- apply(res_nr_db, 2, var)

mean_if <- colMeans(res_if)
mean_if_db <- colMeans(res_if_db)
var_if_db <- apply(res_if_db, 2, var)

mean_ols <- mean(res_ols)
var_ols <- var(res_ols)

mean_pred_p <- colMeans(res_p)
mean_pred_b <- colMeans(res_b)

pdf('2_p_250_Naive-vs-Debiased.pdf', height = 4, width = 12)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2),    # inner margins
    oma = c(0.5, 0.5, 0.5, 0.5),  # outer margins
    mgp = c(3.4, 1.1, 0))
cex_main <- 1.5
cex_lab  <- 1.5
cex_axis <- 1.35

plot(k_all, mean_int - rho, ylim = c(-0.7, 0.7), ylab = 'Bias', xlab = expression(lambda),
     pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     main = 'Integral-based Estimators (Two Splits)')
points(k_all, mean_int_db - rho, col = 'blue', pch = 19, cex = 0.75)
lines(k_all, sum(alpha * beta) * (1 - gammas) - rho * zetas, col = 'red')
lines(k_all, rep(0, length(k_all)), col = 'red')

plot(k_all, mean_nr - rho, ylim = c(-0.7, 0.7), ylab = 'Bias', xlab = expression(lambda),
     pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     main = 'Newey and Robins Estimators (Two Splits)')
points(k_all, mean_nr_db - rho, col = 'blue', pch = 19, cex = 0.75)
lines(k_all, sum(alpha * beta) * (1 - nus), col = 'red')
lines(k_all, rep(0, length(k_all)), col = 'red')

plot(k_all, mean_if - rho, ylim = c(-0.7, 0.7), ylab = 'Bias', xlab = expression(lambda),
     pch = 19, cex = 0.75, cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     main = 'Doubly Robust Estimators (Two Splits)')
points(k_all, mean_if_db - rho, col = 'blue', pch = 19, cex = 0.75)
lines(k_all, sum(alpha * beta) * kappas + rho * zetas, col = 'red')
lines(k_all, rep(0, length(k_all)), col = 'red')

dev.off()


pdf('2_p_250_Variances.pdf', height = 4, width = 12)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2),    # inner margins
    oma = c(0.5, 0.5, 0.5, 0.5),  # outer margins
    mgp = c(3.4, 1.1, 0))
cex_main <- 1.5
cex_lab  <- 1.5
cex_axis <- 1.35

plot(k_all, var_int_db, ylab = 'Variance', xlab = expression(lambda),
      main = 'Integral-based Estimator (Two Splits)',cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75, ylim = c(0, 1))
abline(v = k_all[which.min(var_int_db)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

plot(k_all, var_nr_db, ylab = 'Variance', xlab = expression(lambda),
        main = 'Newey and Robins Estimator (Two Splits)',cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75, ylim = c(0, 1))
abline(v = k_all[which.min(var_nr_db)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

plot(k_all, var_if_db, ylab = 'Variance', xlab = expression(lambda),
     main = 'Doubly Robust Estimator (Two Splits)',cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.75, ylim = c(0, 1))
abline(v = k_all[which.min(var_if_db)], col = 'red', lty = 2)
abline(v = k_all[which.min(mean_pred_p)], col = 'blue', lty = 2)

dev.off()


pdf('p_250_Pred-k.pdf', height = 5, width = 10)
par(mfrow = c(1, 2))

plot(k_all, mean_pred_p, xlab = expression(lambda),ylab = 'MSE', main = 'Estimator of p', cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.5)
lines(k_all, sum(alpha_sim^2) * (1 + gammas - 2 * gammas_p) + zetas, col = 'red')
abline(v = k_all[which.min(mean_pred_p)], col = 'red', lty = 2)

plot(k_all, mean_pred_b, xlab = expression(lambda),ylab = 'MSE', main = 'Estimator of b', cex.lab = cex_lab, cex.axis = cex_axis, cex.main = cex_main,
     pch = 19, cex = 0.5)
lines(k_all, sum(alpha_sim^2) * (1 + gammas - 2 * gammas_p) + zetas, col = 'red')
abline(v = k_all[which.min(mean_pred_b)], col = 'red', lty = 2)
dev.off()
