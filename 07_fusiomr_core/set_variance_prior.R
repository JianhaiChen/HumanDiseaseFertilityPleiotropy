# Set variance priors for FusioMRs (both sigma^2_gamma and sigma^2_theta)
# LD-free; uses clumped, ~independent SNPs
#
# Inputs:
# ghat, gse   : QTL/eQTL effects and SEs (exposure)
# Ghat, Gse   : GWAS effects and SEs (outcome)
# beta0       : optional initial beta (if NULL, GLS init on outcome SEs)
# K           : number of IVs after filtering
# Kmin, Kmax  : floor and cap of K when calculating a
# rho_ov      : sampling correlation due to overlap in [-1,1] (default 0)
# c_gamma     : prior weight per IV for sigma^2_gamma (a_gamma = 1 + c_gamma*K/2)
# c_theta     : prior weight per IV for sigma^2_theta (a_theta = 1 + c_theta*K/2)
# global_mean_gamma, global_mean_theta : global EB centers (for hybrid mode)
# hybrid      : if TRUE, prior mean = eta*local + (1-eta)*global (eta = K/(K+kappa_hybrid))
# kappa_hybrid : pooling control (>0)
# z_thresh    : optional |Z_gamma| selection threshold used to pick QTLs (winner’s-curse fix).
#             Example: for p=5e-8 two-sided, use qnorm(1 - 5e-8/2).
# trim        : tail prob. for winsor
# kappa_gamma : tunning parameter for prior mean of sigma2_gamma 
# kappa_theta : tunning parameter for prior mean of sigma2_theta 
#
# Returns: list with sublists $gamma and $theta:
# $gamma/$theta: list(mom_local, a, b, prior_mean)
# plus shared K, eta, beta0, notes

set_variance_priors <- function(ghat, gse, Ghat, Gse, beta0 = NULL, K = NULL, Kmin = 5, Kmax = 20, 
rho_ov = 0, c_gamma = 0.5, c_theta = 0.8, global_mean_gamma = NULL, global_mean_theta = NULL, 
hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1, kappa_gamma = 1, kappa_theta = 1
){
stopifnot(length(ghat) == length(gse), length(Ghat) == length(Gse), length(ghat) == length(Ghat))
if (is.null(K)) K = length(ghat) # summary statistics are after filtering
notes = c()
# helper: truncated-normal inflation for Z^2 under two-sided |Z|>c selection
kappa_TN <- function(c) {
if (is.null(c) || c <= 0) return(1)
tail = 1 - pnorm(c)
1 + c * dnorm(c) / tail
}
infl_g = kappa_TN(z_thresh)

winsor <- function(x, p = trim) {
if (length(x) < 5) return(x)
qs = quantile(x, c(p, 1-p), na.rm = TRUE)
pmin(pmax(x, qs[[1]]), qs[[2]])
}
wmean <- function(x, p = trim) mean(winsor(x, p), na.rm = TRUE)
wmed  <- function(x, p = trim) median(winsor(x, p), na.rm = TRUE)

if (!is.null(z_thresh) && z_thresh > 0)
notes = c(notes, sprintf("TN correction for gamma with |Z|>%.2f (infl=%.3f)", z_thresh, infl_g))

# sigma^2_gamma: local MoM (de-noised, TN-corrected if requested)
# gamma: winsorized mean(z^2) - 1, rescaled by robust se^2 scale
zg2 = (ghat / gse)^2
mom_gamma = wmed(gse^2, trim) * max(wmean(zg2, trim) - 1, 0)
#Zg = ghat / gse; mom_gamma = mean((gse * Zg)^2 - (gse^2) * infl_g); mom_gamma = max(mom_gamma, 0)

# beta0 (for theta); GLS init if not provided
if (is.null(beta0)) {
wG = 1 / (Gse^2)
xtx = sum((ghat^2) * wG); xty = sum((ghat * Ghat) * wG)
beta0 = if (xtx > 0) xty/xtx else 0
notes = c(notes, "beta0 estimated by GLS on outcome SEs")
} else {
notes = c(notes, "beta0 provided by user")
}
# sigma^2_theta: residual MoM (de-noised, overlap-aware, TN-correct gamma part)
r = Ghat - beta0 * ghat
sv = (Gse^2) + (beta0^2) * (gse^2) - 2 * beta0 * rho_ov * (Gse * gse)
if (!is.null(z_thresh) && z_thresh > 0) {
# conservative: apply TN correction to the gamma-driven sampling component
sv = sv + (beta0^2) * gse^2 * (infl_g - 1)
notes = c(notes, "TN correction applied to gamma component inside residual variance")
}
zr2 = (r^2) / sv
mom_theta = wmed(sv, trim) * max(wmean(zr2, trim) - 1, 0)

#ZG = Ghat / Gse; r = (Gse * ZG) - beta0 * (gse * Zg)
#samp_var_r = (Gse^2) + (beta0^2) * (gse^2) - 2 * beta0 * rho_ov * (Gse * gse)
#if (!is.null(z_thresh) && z_thresh > 0) {
## conservative: apply TN correction to the gamma-driven sampling component
#samp_var_r = samp_var_r + (beta0^2) * gse^2 * (infl_g - 1)
#notes = c(notes, "TN correction applied to gamma component inside residual variance")
#}
#mom_theta = mean( r^2 - samp_var_r )
#mom_theta = max(mom_theta, 0)
  
# shapes (pseudo-df) and hybrid prior means
a_gamma = 1 + c_gamma * min(max(K, Kmin), Kmax) / 2
a_theta = 1 + c_theta * min(max(K, Kmin), Kmax) / 2
if (a_gamma <= 1) a_gamma = 1 + 1e-3
if (a_theta <= 1) a_theta = 1 + 1e-3

if (isTRUE(hybrid)) {
if (is.null(global_mean_gamma) || is.null(global_mean_theta))
stop("hybrid=TRUE requires global_mean_gamma and global_mean_theta.")
eta = K / (K + kappa_hybrid)
m_gamma = eta * mom_gamma + (1 - eta) * global_mean_gamma
m_theta = eta * mom_theta + (1 - eta) * global_mean_theta
notes = c(notes, sprintf("hybrid prior means with eta=K/(K+kappa)=%.3f", eta))
} else {
eta = NA
m_gamma = mom_gamma
m_theta = mom_theta
}

# conservative inflation 
m_gamma = m_gamma * kappa_gamma
m_theta = m_theta * kappa_theta

# scales b to make E[sigma^2] = prior mean (valid for a>1)
b_gamma = (a_gamma - 1) * m_gamma
b_theta = (a_theta - 1) * m_theta
  
# assemble outputs
list(
gamma = list(mom_local  = mom_gamma, a = a_gamma, b = b_gamma, prior_mean = b_gamma / (a_gamma - 1)),
theta = list(mom_local  = mom_theta, a = a_theta, b = b_theta, prior_mean = b_theta / (a_theta - 1)),
K = K, eta = eta, beta0 = beta0, notes  = unique(notes)
)
}


#==============================================================================================


# Set variance priors for FusioMRm (shared exposure, 2 outcomes), LD-free
#
# Inputs:
# ghat, gse          : vectors (length K) of QTL effects and SEs (exposure)
# Ghat_mat, Gse_mat  : K x 2 matrices of GWAS effects and SEs for two outcomes (cols 1,2)
# beta0              : optional length-2 vector of initial betas (if NULL, GLS init per outcome)
# K, Kmin, Kmax      : number of IVs after filtering (auto from ghat if NULL); floor/cap for a/nu
# rho12              : sampling correlation between outcome 1 and 2 (due to sample overlap)
# rho1g, rho2g       : sampling correlations between outcome j and QTL (default 0)
# c_gamma, c_theta   : prior weight per IV for sigma^2_gamma and Sigma_theta
# global_mean_gamma  : optional scalar global EB center for sigma^2_gamma (for hybrid)
# global_mean_theta  : optional 2x2 global EB center for Sigma_theta (for hybrid)
# hybrid             : if TRUE, prior mean = eta*local + (1-eta)*global, eta = K/(K + kappa_hybrid)
# kappa_hybrid       : pooling control (>0)
# z_thresh           : optional |Z_gamma| threshold used to select QTLs (TN winner’s-curse fix)
# trim               : tail prob. for robust trimmed mean/median (set 0 for no trimming)
# kappa_gamma        : tunning parameter for prior mean of sigma2_gamma 
# kappa_theta        : tunning parameter for prior mean of sigma2_theta 
#
# Returns:
# list(
# gamma = list(mom_local, a, b, prior_mean),
# theta = list(mom_local (2x2), nu, Phi, prior_mean (2x2)),
# K, eta, beta0 (length-2), notes
# )
set_variance_priors_m2 <- function(ghat, gse, Ghat_mat, Gse_mat, beta0 = NULL, K = NULL, 
Kmin = 5, Kmax = 20, rho12 = 0, rho1g = 0, rho2g = 0, c_gamma = 0.5, c_theta = 0.8, 
global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, 
z_thresh = NULL, trim = 0.10, kappa_gamma = 1, kappa_theta = 1
) {
  
stopifnot(is.matrix(Ghat_mat), is.matrix(Gse_mat), ncol(Ghat_mat) == 2, ncol(Gse_mat) == 2)
stopifnot(length(ghat) == length(gse), nrow(Ghat_mat) == length(ghat), nrow(Gse_mat) == length(ghat))
if (is.null(K)) K = length(ghat)
notes = c()
  
# helpers 
# TN inflation E[Z^2 | |Z|>c]
kappa_TN <- function(c) {              
if (is.null(c) || c <= 0) return(1)
tail = 1 - pnorm(c)
1 + c * dnorm(c) / tail
}
# winsorize
winsor <- function(x, p) {
if (p <= 0 || length(x) < 5) return(x)
qs = stats::quantile(x, c(p, 1 - p), na.rm = TRUE)
pmin(pmax(x, qs[[1]]), qs[[2]])
}
wmean <- function(x, p) mean(winsor(x, p), na.rm = TRUE)
wmed  <- function(x, p) stats::median(winsor(x, p), na.rm = TRUE)
# PSD projection (eigen floor at 0)
proj_psd <- function(M) {              
M = as.matrix(M); M = (M + t(M)) / 2
ev = eigen(M, symmetric = TRUE)
vals = pmax(ev$values, 0)
ev$vectors %*% diag(vals, nrow = 2) %*% t(ev$vectors)
}
  
infl_g = kappa_TN(z_thresh)
if (!is.null(z_thresh) && z_thresh > 0)
notes = c(notes, sprintf("TN correction for gamma with |Z|>%.2f (infl=%.3f)", z_thresh, infl_g))

# sigma^2_gamma, local MoM (scalar, de-noised, TN-corrected if requested)
Zg = ghat / gse
# de-noised MoM; robust version (trim) or exact subtraction with TN factor
if (trim > 0) {
mom_gamma = wmed(gse^2, trim) * max(wmean(Zg^2, trim) - 1, 0)
} else {
mom_gamma = mean((gse * Zg)^2 - gse^2 * infl_g)
mom_gamma = max(mom_gamma, 0)
}
  
# beta0 for each outcome (GLS init if missing)
if (is.null(beta0)) {
beta0 = numeric(2)
for (j in 1:2) {
wj = 1 / (Gse_mat[, j]^2)
keep = which(is.finite(wj))
xtx = sum(ghat[keep]^2 * wj[keep])
xty = sum(ghat[keep] * Ghat_mat[, j][keep] * wj[keep])
beta0[j] = if (xtx > 0) xty / xtx else 0
}
notes = c(notes, "beta0 estimated by GLS per outcome")
} else {
stopifnot(length(beta0) == 2)
notes = c(notes, "beta0 provided by user")
}
  
# Sigma_theta (2x2) local MoM 
# residuals r_k (2x1), sampling covariance S_k (2x2) with overlap and TN fix on gamma parts
sg = gse
s1 = Gse_mat[, 1]; s2 = Gse_mat[, 2]
b1 = beta0[1]; b2 = beta0[2]
Z1 = Ghat_mat[, 1] / s1; Z2 = Ghat_mat[, 2] / s2
r1 = (s1*Z1) - b1 * (sg*Zg)
r2 = (s2*Z2) - b2 * (sg*Zg)
  
S11 = (s1^2) + (b1^2) * (sg^2) - 2 * b1 * rho1g * (s1*sg)
S22 = (s2^2) + (b2^2) * (sg^2) - 2 * b2 * rho2g * (s2*sg)
S12 = rho12 * (s1*s2) + b1*b2 * (sg^2) - b1 * rho2g * (s2*sg) - b2 * rho1g * (s1*sg)
if (!is.null(z_thresh) && z_thresh > 0) {
S11 = S11 + (b1^2) * sg^2 * (infl_g-1)
S22 = S22 + (b2^2) * sg^2 * (infl_g-1)
S12 = S12 + (b1*b2)* sg^2 * (infl_g-1)
notes = c(notes, "TN correction applied to gamma components inside residual sampling covariance")
}
  
# local MoM covariance: average of (r r^T - S_k), PSD-projected
M11 = mean(r1*r1 - S11)
M22 = mean(r2*r2 - S22)
M12 = mean(r1*r2 - S12)
mom_theta_mat = matrix(c(M11, M12, M12, M22), 2, 2)
mom_theta_mat = proj_psd(mom_theta_mat)

# pseudo-df and hybrid centers 
Keff = min(max(K, Kmin), Kmax)
a_gamma = 1 + (c_gamma * Keff) / 2
# For IW(p=2), need nu > p+1 = 3 for mean to exist
nu_theta = 3 + (c_theta * Keff)
if (a_gamma <= 1) a_gamma = 1 + 1e-3
if (nu_theta <= 3) nu_theta = 3 + 1e-3
  
if (isTRUE(hybrid)) {
if (is.null(global_mean_gamma) || is.null(global_mean_theta))
stop("hybrid=TRUE requires global_mean_gamma (scalar) and global_mean_theta (2x2).")
if (!is.matrix(global_mean_theta) || any(dim(global_mean_theta) != c(2, 2)))
stop("global_mean_theta must be a 2x2 matrix.")
eta = K / (K + kappa_hybrid)
m_gamma = eta * mom_gamma + (1 - eta) * global_mean_gamma
m_theta = eta * mom_theta_mat + (1 - eta) * global_mean_theta
notes = c(notes, sprintf("hybrid prior means with eta=K/(K+kappa)=%.3f", eta))
} else {
eta = NA_real_
m_gamma = mom_gamma
m_theta = mom_theta_mat
}

# conservative inflation 
m_gamma = m_gamma * kappa_gamma
m_theta = m_theta * kappa_theta 

# ensure m_theta is PSD
m_theta = proj_psd(m_theta)

# scales so that E[sigma^2_gamma] = m_gamma and E[Sigma_theta] = m_theta
b_gamma = (a_gamma - 1) * m_gamma
Phi = (nu_theta - 3) * m_theta
  
# outputs 
list(
gamma = list(mom_local = mom_gamma, a = a_gamma, b = b_gamma, prior_mean = b_gamma / (a_gamma - 1)),
theta = list(mom_local = mom_theta_mat, nu = nu_theta, Phi = Phi, prior_mean = Phi / (nu_theta - 3)),
K = K, eta = eta, beta0 = beta0, notes = unique(notes)
)
}


#==============================================================================================



# Set variance priors for FusioMRm (two exposures, two outcomes), diagonal B, LD-free
#
# Inputs:
# ghat_mat, gse_mat : K x 2 effects and SEs for two exposures 
# Ghat_mat, Gse_mat : K x 2 effects and SEs for two outcomes 
# B0                : optional length-2 initial betas (beta1, beta2). If NULL, GLS init per outcome.
# K, Kmin, Kmax     : #IVs (auto from nrow if NULL); floor/cap for prior strength
# Overlap correlations (sampling):
# rho12             : outcome1–outcome2 overlap (default 0)
# rho_gg            : exposure1–exposure2 QTL overlap (default 0)
# rho_gj            : list length 2 with exposure–outcome overlaps:
#                     rho_gj[[1]] = c(rho_{γ1,1}, rho_{γ1,2})
#                     rho_gj[[2]] = c(rho_{γ2,1}, rho_{γ2,2})
# Hyper-parameters:
# c_gamma, c_theta  : prior weight per IV for IW on Sigma_gamma, Sigma_theta
# global_mean_gamma : optional 2x2 global EB center for Sigma_gamma (for hybrid)
# global_mean_theta : optional 2x2 global EB center for Sigma_theta (for hybrid)
# hybrid            : if TRUE, prior mean = eta*local + (1-eta)*global, eta = K/(K + kappa_hybrid)
# kappa_hybrid      : pooling control (>0)
# Robust/selection knobs:
# z_thresh          : optional |Z_gamma| threshold (TN winner’s-curse fix)
# Conservative inflation:
# kappa_gamma       : scalar multiplier for prior mean of Sigma_gamma (>=1)
# kappa_theta       : scalar multiplier for prior mean of Sigma_theta (>=1)
#
# Returns:
# list(
# gamma = list(mom_local (2x2), nu, Phi, prior_mean (2x2)),
# theta = list(mom_local (2x2), nu, Phi, prior_mean (2x2)),
# K, Keff, eta, B0 (length-2), notes
# )

set_variance_priors_m2x2_diag <- function(ghat_mat, gse_mat, Ghat_mat, Gse_mat, B0 = NULL, K = NULL, 
Kmin = 5, Kmax = 20, rho12 = 0, rho_gg = 0, rho_gj = list(c(0,0), c(0,0)), c_gamma = 0.8, c_theta = 0.8,
global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL,
kappa_gamma = 1, kappa_theta = 1) 
{
## --- checks ---
stopifnot(is.matrix(ghat_mat), is.matrix(gse_mat), is.matrix(Ghat_mat), is.matrix(Gse_mat))
stopifnot(ncol(ghat_mat) == 2, ncol(gse_mat) == 2, ncol(Ghat_mat) == 2, ncol(Gse_mat) == 2)
stopifnot(nrow(ghat_mat) == nrow(gse_mat), nrow(Ghat_mat) == nrow(ghat_mat), nrow(Gse_mat) == nrow(ghat_mat))
if (is.null(K)) K = nrow(ghat_mat)
notes = c()

## --- helpers ---
kappa_TN <- function(c) { # E[Z^2 | |Z|>c] for Z~N(0,1)
if (is.null(c) || c <= 0) return(1)
tail = 1 - pnorm(c)
1 + c * dnorm(c) / tail
}

proj_psd <- function(M) { # PSD projection (eigen floor at 0)
M = (M + t(M))/2
ev = eigen(M, symmetric = TRUE)
V = ev$vectors; d <- pmax(ev$values, 0)
V %*% diag(d, nrow = 2) %*% t(V)
}

make_spd <- function(M, eps = 1e-8) {
M = (M + t(M))/2
e = eigen(M, symmetric = TRUE)
d = pmax(e$values, eps * mean(pmax(e$values, 0) + 1e-12))
e$vectors %*% diag(d, 2) %*% t(e$vectors)
}

infl_g = kappa_TN(z_thresh)
if (!is.null(z_thresh) && z_thresh > 0)
notes = c(notes, sprintf("TN correction for gamma with |Z|>%.2f (infl=%.3f)", z_thresh, infl_g))

## --- local MoM for Sigma_gamma (2x2 across exposures) ---
s_g1 = gse_mat[,1]; s_g2 = gse_mat[,2]
Sg11 = s_g1^2; Sg22 = s_g2^2; Sg12 = rho_gg * (s_g1 * s_g2)
if (!is.null(z_thresh) && z_thresh > 0) {
Sg11 = Sg11 * infl_g; Sg22 = Sg22 * infl_g; Sg12 = Sg12 * infl_g
notes = c(notes, "TN correction applied inside S^(gamma)")
}
G1 = ghat_mat[,1]; G2 = ghat_mat[,2]
Mg11 = mean(G1*G1 - Sg11)
Mg22 = mean(G2*G2 - Sg22)
Mg12 = mean(G1*G2 - Sg12)
mom_gamma_mat = matrix(c(Mg11, Mg12, Mg12, Mg22), 2, 2)
mom_gamma_mat = proj_psd(mom_gamma_mat)
mom_gamma_mat = kappa_gamma * mom_gamma_mat # optional conservative inflation

## --- initial GLS for diagonal B0 if missing ---
if (is.null(B0)) {
B0 = numeric(2)
ridge = 1e-12
# outcome 1 on exposure 1
w1 = 1 / (Gse_mat[,1]^2)
xtx1 = sum(ghat_mat[,1]^2 * w1) + ridge
xty1 = sum(ghat_mat[,1] * Ghat_mat[,1] * w1)
B0[1] = if (xtx1 > 0) xty1 / xtx1 else 0
# outcome 2 on exposure 2
w2 = 1 / (Gse_mat[,2]^2)
xtx2 = sum(ghat_mat[,2]^2 * w2) + ridge
xty2 = sum(ghat_mat[,2] * Ghat_mat[,2] * w2)
B0[2] = if (xtx2 > 0) xty2 / xtx2 else 0
notes = c(notes, "B0 (beta1,beta2) estimated by GLS per outcome on its matching exposure")
} else {
stopifnot(is.numeric(B0), length(B0) == 2)
notes = c(notes, "B0 provided by user (diagonal model)")
}

## --- local MoM for Sigma_theta (2x2 across outcomes), diagonal B ---
# residuals r_k = Γ_hat - diag(B0) * γ_hat (matching exposure only)
R1 = Ghat_mat[,1] - (B0[1] * ghat_mat[,1])
R2 = Ghat_mat[,2] - (B0[2] * ghat_mat[,2])

# sampling covariance S^(θ)_k with overlaps
s1 = Gse_mat[,1]; s2 = Gse_mat[,2]
r_g1_1 = rho_gj[[1]][1]; r_g1_2 = rho_gj[[1]][2]  # γ1 with outcome1, outcome2
r_g2_1 = rho_gj[[2]][1]; r_g2_2 = rho_gj[[2]][2]  # γ2 with outcome1, outcome2
S11 = (s1^2) + (B0[1]^2) * (s_g1^2) - 2 * B0[1] * r_g1_1 * (s_g1 * s1)
S22 = (s2^2) + (B0[2]^2) * (s_g2^2) - 2 * B0[2] * r_g2_2 * (s_g2 * s2)
# Covariance term includes outcome-outcome overlap, two cross exposure–outcome overlaps, and QTL overlap
S12 = rho12 * (s1 * s2) - B0[1] * r_g1_2 * (s_g1 * s2) - B0[2] * r_g2_1 * (s_g2 * s1) + (B0[1] * B0[2]) * (rho_gg * s_g1 * s_g2)

# local MoM covariance: average of (r r^T - S_k), then PSD-project
M11 = mean(R1*R1 - S11)
M22 = mean(R2*R2 - S22)
M12 = mean(R1*R2 - S12)
mom_theta_mat = matrix(c(M11, M12, M12, M22), 2, 2)
mom_theta_mat = proj_psd(mom_theta_mat)
mom_theta_mat = kappa_theta * mom_theta_mat # optional conservative inflation

## --- pseudo-df and hybrid centers ---
Keff <- min(max(K, Kmin), Kmax)
# IW(p=2): need nu > p+1 = 3
nu_gamma = 3 + (c_gamma * Keff)
nu_theta = 3 + (c_theta * Keff)
if (nu_gamma <= 3) nu_gamma = 3 + 1
if (nu_theta <= 3) nu_theta = 3 + 1

if (isTRUE(hybrid)) {
if (is.null(global_mean_gamma) || is.null(global_mean_theta))
stop("hybrid=TRUE requires global_mean_gamma (2x2) and global_mean_theta (2x2).")
stopifnot(is.matrix(global_mean_gamma), all(dim(global_mean_gamma) == c(2,2)))
stopifnot(is.matrix(global_mean_theta), all(dim(global_mean_theta) == c(2,2)))
eta = K / (K + kappa_hybrid)
m_gamma = eta * mom_gamma_mat + (1 - eta) * global_mean_gamma
m_theta = eta * mom_theta_mat + (1 - eta) * global_mean_theta
notes = c(notes, sprintf("hybrid prior means with eta=K/(K+kappa)=%.3f (diagonal B)", eta))
} else {
eta = NA_real_
m_gamma = mom_gamma_mat
m_theta = mom_theta_mat
}
m_gamma = proj_psd(m_gamma)
m_theta = proj_psd(m_theta)
# upgrade PSD -> SPD for correct inverse
m_gamma = make_spd(m_gamma)   
m_theta = make_spd(m_theta)

## --- set IW scales so that E[Sigma_*] = prior mean ---
Phi_gamma = (nu_gamma - 3) * m_gamma
Phi_theta = (nu_theta - 3) * m_theta

list(
gamma = list(mom_local  = mom_gamma_mat, nu = nu_gamma, Phi = Phi_gamma, prior_mean = Phi_gamma / (nu_gamma - 3)),
theta = list(mom_local  = mom_theta_mat, nu = nu_theta, Phi = Phi_theta, prior_mean = Phi_theta / (nu_theta - 3)),
K = K, Keff = Keff, eta = eta, B0 = B0, notes = unique(notes)
)
}


#==============================================================================================

# df: data frame including all gene-SNP pairs (rows) before p-filtering

get_global_mom <- function(df, q_gamma = 0.5, q_theta = 0.5, lb = 1e-6, pos = TRUE) {
# minimal column checks
req = c("gene_symbol", "beta", "se", "b_out", "se_out")
miss = setdiff(req, colnames(df))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
by_gene = df %>% group_by(gene_symbol) %>% group_split()
per_gene = lapply(by_gene, function(d) {
ghat = d$beta; gse = d$se; Ghat = d$b_out; Gse = d$se_out
# sigma^2_gamma: local MoM (de-noised)
Zg = ghat / gse
mom_gamma = mean((gse * Zg)^2 - (gse^2))
mom_gamma = max(mom_gamma, 0)
# beta0 (for theta); GLS init 
wG = 1 / (Gse^2)
beta0 = as.numeric( sum(ghat * Ghat * wG) / sum(ghat^2 * wG) )
# sigma^2_theta: residual MoM (de-noised)
ZG = Ghat / Gse
r = (Gse * ZG) - beta0 * (gse * Zg)
samp_var_r = (Gse^2) + (beta0^2) * (gse^2)
mom_theta = mean(r^2 - samp_var_r)
mom_theta = max(mom_theta, 0)
c(mom_gamma, mom_theta)
})
M = do.call(rbind, per_gene)

safe_quantile <- function(x, p, lb, pos = TRUE) {
if (all(is.na(x))) return(NA_real_)
if (pos) {
qq = as.numeric(quantile(x[x>0], probs = p, na.rm = TRUE, names = FALSE))
} else {
qq = as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}
max(qq, lb)
}

c(
mom_gamma_q = safe_quantile(M[, 1], q_gamma, lb, pos), 
mom_theta_q = safe_quantile(M[, 2], q_theta, lb, pos)
)
}




