
# function to set up initial values for FusioMR method
# se-so, uhp only

init_setup = function(niter, K, alpha_init, beta_init, sigma_gamma_init, sigma_theta_init) {
# starting values
theta_init = rep(0, K) 
gamma_init = rep(0, K)
# MCMC tracks
theta_tk = gamma_tk = matrix(NA, nrow = niter, ncol = K)
alpha_tk = beta_tk = sigma2_gamma_tk = sigma2_theta_tk = rep(NA, niter)
theta_tk[1, ] = theta_init
gamma_tk[1, ] = gamma_init
alpha_tk[1] = alpha_init
beta_tk[1] = beta_init
sigma2_gamma_tk[1] = sigma_gamma_init^2
sigma2_theta_tk[1] = sigma_theta_init^2
# output
res = list(theta_tk = theta_tk, gamma_tk = gamma_tk, 
alpha_tk = alpha_tk, beta_tk = beta_tk,  
sigma2_gamma_tk = sigma2_gamma_tk, sigma2_theta_tk = sigma2_theta_tk)
return(res)
}


