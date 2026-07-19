
# Single-population MR using hierarchical Bayesian model 
# UHP only, no CHP effect

library(invgamma)

gibbs_single_uhp_only = function(niter, K, beta_tk, gamma_tk, theta_tk, 
sigma2_gamma_tk, sigma2_theta_tk, Gamma_hat, gamma_hat, s2_hat_Gamma, s2_hat_gamma, 
a_gamma, b_gamma, a_theta, b_theta) {
  
# Set current values as initial values
beta_cur = beta_tk[1]
theta_cur = theta_tk[1, ]
gamma_cur = gamma_tk[1, ]
sigma2_gamma_cur = sigma2_gamma_tk[1]
sigma2_theta_cur = sigma2_theta_tk[1]

for (iter in 1:(niter - 1)) {

# Update gamma_k
for (k in 1:K) {
Ak_gamma = beta_cur^2 / s2_hat_Gamma[k] + 1 / s2_hat_gamma[k] + 1 / sigma2_gamma_cur
Bk_gamma = beta_cur * (Gamma_hat[k] - theta_cur[k]) / s2_hat_Gamma[k] + gamma_hat[k] / s2_hat_gamma[k]
gamma_cur[k] = rnorm(n = 1, mean = Bk_gamma / Ak_gamma, sd = sqrt(1 / Ak_gamma))
}
gamma_tk[iter + 1, ] = gamma_cur
  
# Update theta_k
for (k in 1:K) {
Ak_theta = 1 / s2_hat_Gamma[k] + 1 / sigma2_theta_cur
Bk_theta = (Gamma_hat[k] - beta_cur * gamma_cur[k]) / s2_hat_Gamma[k]
theta_cur[k] = rnorm(n = 1, mean = Bk_theta / Ak_theta, sd = sqrt(1 / Ak_theta))
}
theta_tk[iter + 1, ] = theta_cur

# Update sigma2_gamma and sigma2_theta
sigma2_gamma_cur = invgamma::rinvgamma(n = 1, shape = a_gamma + K/2, rate = b_gamma + 0.5 * sum(gamma_cur^2))
sigma2_theta_cur = invgamma::rinvgamma(n = 1, shape = a_theta + K/2, rate = b_theta + 0.5 * sum(theta_cur^2))
if (sigma2_gamma_cur < 1e-10) sigma2_gamma_cur = 1e-10
if (sigma2_gamma_cur > 10) sigma2_gamma_cur = 10
if (sigma2_theta_cur < 1e-10) sigma2_theta_cur = 1e-10
if (sigma2_theta_cur > 10) sigma2_theta_cur = 10
sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur
sigma2_theta_tk[iter + 1] = sigma2_theta_cur

# Update beta
U_beta = Gamma_hat - theta_cur 
Omega_hat_Gamma = diag(1/s2_hat_Gamma)
# A_beta = t(gamma_cur) %*% Omega_hat_Gamma %*% gamma_cur
A_beta = max(t(gamma_cur) %*% Omega_hat_Gamma %*% gamma_cur, 1e-6)
mu_beta = t(gamma_cur) %*% Omega_hat_Gamma %*% U_beta / A_beta
beta_cur = rnorm(n = 1, mean = mu_beta, sd = sqrt(1 / A_beta))
beta_tk[iter + 1] = beta_cur
} # end iter

# output
out = list(K = K , beta_tk = beta_tk, gamma_tk = gamma_tk, theta_tk = theta_tk, 
sigma2_gamma_tk = sigma2_gamma_tk, sigma2_theta_tk = sigma2_theta_tk)
return(out)
}

