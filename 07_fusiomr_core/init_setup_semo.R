
# function to set up initial values for FusioMR method
# se-mo, uhp only

init_setup_semo_uhp_only = function(niter, K, beta_1_init, beta_2_init, sigma_gamma_init) {
# starting values
theta_1_init = theta_2_init = rep(0, K) 
gamma_init = rep(0, K)
# MCMC tracks
gamma_tk = theta_1_tk = theta_2_tk = matrix(NA, nrow = niter, ncol = K)
beta_1_tk = beta_2_tk = sigma2_gamma_tk = rep(NA, niter)
gamma_tk[1, ] = gamma_init
theta_1_tk[1, ] = theta_1_init
theta_2_tk[1, ] = theta_2_init
beta_1_tk[1] = beta_1_init
beta_2_tk[1] = beta_2_init
sigma2_gamma_tk[1] = sigma_gamma_init^2
# output
res = list(theta_1_tk = theta_1_tk, theta_2_tk = theta_2_tk, gamma_tk = gamma_tk, 
beta_1_tk = beta_1_tk, beta_2_tk = beta_2_tk,
sigma2_gamma_tk = sigma2_gamma_tk)
return(res)
}


