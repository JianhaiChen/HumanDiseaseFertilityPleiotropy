#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// rinvgamma from R::invgamma
// [[Rcpp::export]]
NumericVector my_rinvgamma(int n, double shape, double rate) {
Function ff("rinvgamma", Environment::namespace_env("invgamma"));
NumericVector res = ff(n, Named("shape") = shape, _["rate"] = rate);
return res;
}

// [[Rcpp::export]]
List gibbs_seso_uhp_only_cpp(int niter, int K, 
NumericVector beta_tk, NumericMatrix theta_tk, NumericMatrix gamma_tk, 
NumericVector sigma2_gamma_tk, NumericVector sigma2_theta_tk,
NumericVector Gamma_hat, NumericVector gamma_hat, 
NumericVector s2_hat_Gamma, NumericVector s2_hat_gamma, 
double a_gamma, double b_gamma, double a_theta, double b_theta) 
{
// Initialize values
double beta_cur = beta_tk[0];
NumericVector theta_cur = theta_tk(0, _);
NumericVector gamma_cur = gamma_tk(0, _);
double sigma2_gamma_cur = sigma2_gamma_tk[0];
double sigma2_theta_cur = sigma2_theta_tk[0];

// iteration
for (int iter = 0; iter < (niter-1); iter ++) {
// Update gamma_k
for (int k = 0; k < K; k++) {
double Ak_gamma = beta_cur * beta_cur / s2_hat_Gamma[k] + 1.0 / s2_hat_gamma[k] + 1.0 / sigma2_gamma_cur;
double Bk_gamma = beta_cur * (Gamma_hat[k] - theta_cur[k]) / s2_hat_Gamma[k] + gamma_hat[k] / s2_hat_gamma[k];
gamma_cur[k] = rnorm(1, Bk_gamma / Ak_gamma, sqrt(1.0 / Ak_gamma))[0];
}
gamma_tk(iter + 1, _) = gamma_cur;
// Update theta_k
for (int k = 0; k < K; k++) {
double Ak_theta = 1.0 / s2_hat_Gamma[k] + 1.0 / sigma2_theta_cur;
double Bk_theta = (Gamma_hat[k] - beta_cur * gamma_cur[k]) / s2_hat_Gamma[k];
theta_cur[k] = rnorm(1, Bk_theta / Ak_theta, sqrt(1.0 / Ak_theta))[0];
}
theta_tk(iter + 1, _) = theta_cur;
// Update sigma2_gamma and sigma2_theta
sigma2_gamma_cur = my_rinvgamma(1, a_gamma + K / 2.0, b_gamma + 0.5 * sum(pow(gamma_cur, 2)))[0];
sigma2_theta_cur = my_rinvgamma(1, a_theta + K / 2.0, b_theta + 0.5 * sum(pow(theta_cur, 2)))[0];
//if (sigma2_gamma_cur < 1e-10) sigma2_gamma_cur = 1e-10;
//if (sigma2_gamma_cur > 10) sigma2_gamma_cur = 10;
//if (sigma2_theta_cur < 1e-10) sigma2_theta_cur = 1e-10;
//if (sigma2_theta_cur > 10) sigma2_theta_cur = 10;
sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;
sigma2_theta_tk[iter + 1] = sigma2_theta_cur;
// Update beta
vec U_beta = Gamma_hat - theta_cur;
vec W_beta = gamma_cur;
mat Omega_hat_Gamma(K, K);  
for (int ii = 0; ii < K; ii ++) {
Omega_hat_Gamma(ii, ii) = 1/s2_hat_Gamma[ii];       
}
mat A_beta  = trans(W_beta) * Omega_hat_Gamma * W_beta;
mat mu_beta = trans(W_beta) * Omega_hat_Gamma * U_beta / A_beta; 
beta_cur = rnorm(1, mu_beta(0,0), sqrt(1/A_beta(0,0)))[0];
beta_tk[iter + 1] = beta_cur;
} // end iter
List res = List::create(
Named("K") = K, Named("beta_tk") = beta_tk, Named("gamma_tk") = gamma_tk, Named("theta_tk") = theta_tk,
Named("sigma2_gamma_tk") = sigma2_gamma_tk, Named("sigma2_theta_tk") = sigma2_theta_tk);
return res;
}


