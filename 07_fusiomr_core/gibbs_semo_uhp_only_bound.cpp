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
NumericVector my_rinvwishart(double nu, mat S) {
Function ff("rinvwishart"); 
NumericVector res = ff(Named("nu") = nu, _["S"] = S);
return res;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
int ncols = sigma.n_cols;
arma::mat Y = arma::randn(n, ncols);
return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// Gibbs sampling
// [[Rcpp::export]]
List gibbs_semo_uhp_only_rcpp(int niter, int K,
vec beta_1_tk, vec beta_2_tk, mat theta_1_tk, mat theta_2_tk, NumericMatrix gamma_tk, NumericVector sigma2_gamma_tk,
vec Gamma_hat_1, vec Gamma_hat_2, vec s2_hat_Gamma_1, vec s2_hat_Gamma_2, vec gamma_hat, vec s2_hat_gamma,
double a_gamma, double b_gamma, mat Sigma_theta_init, double m_theta, mat V_theta) 
{
// Initialize values
vec beta_cur(2); beta_cur[0] = beta_1_tk[0]; beta_cur[1] = beta_2_tk[0];
NumericVector gamma_cur = gamma_tk(0, _);
mat theta_cur(K,2);
for (int ii = 0; ii < K; ii ++) { theta_cur(ii,0) = theta_1_tk(0,ii); theta_cur(ii,1) = theta_2_tk(0,ii); }
double sigma2_gamma_cur = sigma2_gamma_tk[0];
mat Sigma_theta_cur = Sigma_theta_init;
std::vector<double> sigma2_theta1_tk(niter, 0.0);
std::vector<double> sigma2_theta2_tk(niter, 0.0);
sigma2_theta1_tk[0] = Sigma_theta_init(0,0);
sigma2_theta2_tk[0] = Sigma_theta_init(0,0);
  
// iteration
for (int iter = 0; iter < (niter-1); iter ++) {

// Update gamma_k
for (int k = 0; k < K; k++) {
double Ak_gamma = pow(beta_cur[0],2)/s2_hat_Gamma_1[k] + pow(beta_cur[1],2)/s2_hat_Gamma_2[k] + 1.0/s2_hat_gamma[k] + 1.0/sigma2_gamma_cur;
double Bk_gamma = beta_cur[0]*(Gamma_hat_1[k]-theta_cur(k,0))/s2_hat_Gamma_1[k] + beta_cur[1]*(Gamma_hat_2[k]-theta_cur(k,1))/s2_hat_Gamma_2[k] + gamma_hat[k]/s2_hat_gamma[k];
gamma_cur[k] = rnorm(1, Bk_gamma/Ak_gamma, sqrt(1.0/Ak_gamma))[0];
}
gamma_tk(iter + 1, _) = gamma_cur;

// update theta_k
for (int k = 0; k < K; k ++) {
mat S_hat_Gamma_k(2,2);
S_hat_Gamma_k(0,0) = s2_hat_Gamma_1[k]; S_hat_Gamma_k(1,1) = s2_hat_Gamma_2[k];
mat precision_mat_post_theta = inv(S_hat_Gamma_k) + inv(Sigma_theta_cur);
vec temp = {Gamma_hat_1[k] - beta_cur[0]*gamma_cur[k], Gamma_hat_2[k] - beta_cur[1]*gamma_cur[k]};
vec mu_post_nom_theta = inv(S_hat_Gamma_k) * temp;
vec mu_post_theta = inv(precision_mat_post_theta) * mu_post_nom_theta;
mat theta_sample = mvrnormArma(1, mu_post_theta, inv(precision_mat_post_theta));
theta_cur(k,0) = theta_sample(0,0);
theta_cur(k,1) = theta_sample(0,1);
}

// Update sigma2_gamma
sigma2_gamma_cur = my_rinvgamma(1, a_gamma+K/2.0, b_gamma+0.5*sum(pow(gamma_cur,2)))[0];
if (sigma2_gamma_cur < 1e-4) sigma2_gamma_cur = 1e-4;
if (sigma2_gamma_cur > 10) sigma2_gamma_cur = 10;
sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;

//update Sigma_theta
mat S_theta(2,2);
for (int ii = 0; ii < K; ii ++) {
rowvec rowk = theta_cur.row(ii);
S_theta = S_theta + trans(rowk) * rowk;
}
S_theta = S_theta / K;
double m_post_theta = m_theta + K;
mat V_post_theta = S_theta * K + V_theta;
mat Sigma_theta_sample = my_rinvwishart(m_post_theta, V_post_theta);
double rho_theta = Sigma_theta_sample(1,0)/sqrt(Sigma_theta_sample(0,0))/sqrt(Sigma_theta_sample(3,0));
if (Sigma_theta_sample(0,0) < 1e-8) Sigma_theta_sample(0,0) = 1e-8;
if (Sigma_theta_sample(0,0) > 10) Sigma_theta_sample(0,0) = 10;
if (Sigma_theta_sample(3,0) < 1e-8) Sigma_theta_sample(3,0) = 1e-8;
if (Sigma_theta_sample(3,0) > 10) Sigma_theta_sample(3,0) = 10;
Sigma_theta_cur(0,0) = Sigma_theta_sample(0,0);
Sigma_theta_cur(1,1) = Sigma_theta_sample(3,0);
//Sigma_theta_cur(0,1) = Sigma_theta_sample(0,0)*Sigma_theta_sample(3,0)*rho_theta;
Sigma_theta_cur(0,1) = sqrt(Sigma_theta_cur(0,0))*sqrt(Sigma_theta_cur(1,1))*rho_theta;
Sigma_theta_cur(1,0) = Sigma_theta_cur(0,1);
sigma2_theta1_tk[iter+1] = Sigma_theta_cur(0,0);
sigma2_theta2_tk[iter+1] = Sigma_theta_cur(1,1);

// update beta_j
mat Omega_hat_Gamma_1(K, K);
mat Omega_hat_Gamma_2(K, K);
for (int ii = 0; ii < K; ii ++) {
Omega_hat_Gamma_1(ii, ii) = 1/s2_hat_Gamma_1[ii];
Omega_hat_Gamma_2(ii, ii) = 1/s2_hat_Gamma_2[ii];
}
vec U_beta_1(K);
vec U_beta_2(K);
vec W_beta_1(K);
vec W_beta_2(K);
for (int ii = 0; ii < K; ii ++) {
U_beta_1[ii] = Gamma_hat_1[ii] - theta_cur(ii,0);
U_beta_2[ii] = Gamma_hat_2[ii] - theta_cur(ii,1);
W_beta_1[ii] = gamma_cur[ii];
W_beta_2[ii] = gamma_cur[ii];
}
mat A_beta_mat_1  = trans(W_beta_1) * Omega_hat_Gamma_1 * W_beta_1;
mat A_beta_mat_2  = trans(W_beta_2) * Omega_hat_Gamma_2 * W_beta_2;
mat mu_beta_mat_1  = trans(W_beta_1) * Omega_hat_Gamma_1 * U_beta_1;
mat mu_beta_mat_2  = trans(W_beta_2) * Omega_hat_Gamma_2 * U_beta_2;
double A_beta_1 = std::max(A_beta_mat_1(0,0), 1e-6);
double A_beta_2 = std::max(A_beta_mat_2(0,0), 1e-6);
double mu_beta_1 = mu_beta_mat_1(0,0) / A_beta_1; 
double mu_beta_2 = mu_beta_mat_2(0,0) / A_beta_2; 
beta_cur[0] = rnorm(1, mu_beta_1, sqrt(1/A_beta_1))[0];
beta_cur[1] = rnorm(1, mu_beta_2, sqrt(1/A_beta_2))[0];
// mat A_beta_1  = trans(W_beta_1) * Omega_hat_Gamma_1 * W_beta_1;
// mat A_beta_2  = trans(W_beta_2) * Omega_hat_Gamma_2 * W_beta_2;
// mat mu_beta_1 = trans(W_beta_1) * Omega_hat_Gamma_1 * U_beta_1 / A_beta_1;
// mat mu_beta_2 = trans(W_beta_2) * Omega_hat_Gamma_2 * U_beta_2 / A_beta_2;
// beta_cur[0] = rnorm(1, mu_beta_1(0,0), sqrt(1/A_beta_1(0,0)))[0];
// beta_cur[1] = rnorm(1, mu_beta_2(0,0), sqrt(1/A_beta_2(0,0)))[0];
beta_1_tk[iter+1] = beta_cur[0];
beta_2_tk[iter+1] = beta_cur[1];
} // end iter

List res = List::create(Named("K") = K, Named("beta_1_tk") = beta_1_tk, 
Named("beta_2_tk") = beta_2_tk, Named("sigma2_gamma_tk") = sigma2_gamma_tk, 
Named("sigma2_theta1_tk") = sigma2_theta1_tk, 
Named("sigma2_theta2_tk") = sigma2_theta2_tk);
return res;
}



