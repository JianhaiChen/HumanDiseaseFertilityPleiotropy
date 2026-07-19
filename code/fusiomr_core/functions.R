# clumping function
clump <- function(dat, SNP_col = "eQTL_variant_id", pval_col = "rowmeta", clump_kb = 250, clump_r2 = 0.1, clump_p = 0.999, bfile = "/scratch/t.phs.yihaolu/GWAS_summary_statistics/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01", plink_bin = "plink", pop="EUR") {
  df <- data.frame(rsid = dat[, ..SNP_col], pval = dat[,..pval_col])
  colnames(df) = c("rsid", "pval")
  out <- tryCatch({
    ieugwasr::ld_clump(df, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, bfile=bfile, plink_bin = plink_bin, pop = pop)
  }, silent = TRUE, error = function(x) return(NA)
  )
  if(length(out)==1) {
    return(NA)
  }
  MRdat <- dat[which(unlist(dat[,..SNP_col]) %in% out$rsid),]
  return(MRdat)
}

# function to run MR
run_mr = function(gene, b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior) {
  mr.obj = MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
  IVW_f = MendelianRandomization::mr_ivw(mr.obj, model = 'fixed')
  b_ivw_fixed = IVW_f$Estimate; se_ivw_fixed = IVW_f$StdError; p_ivw_fixed = IVW_f$Pvalue
  Egger = tryCatch({MendelianRandomization::mr_egger(mr.obj)}, error = function(e) {NA})
  if(!is.null(Egger) & class(Egger)=='Egger') { 
    b_egger = Egger$Estimate; se_egger = Egger$StdError.Est; p_egger = Egger$Pvalue.Est
  } else b_egger = se_egger = p_egger = NA
  cml = tryCatch({MendelianRandomization::mr_cML(mr.obj, MA = TRUE, DP = FALSE, n = 100)}, error = function(e) {NA})
  if(!is.null(cml)) {  
    b_cml = cml$Estimate; se_cml = cml$StdError; p_cml = cml$Pvalue
  } else b_cml = se_cml = p_cml = NA
  out_other_methods = c(length(b_exp), b_ivw_fixed, se_ivw_fixed, p_ivw_fixed, b_egger, se_egger, p_egger, b_cml, se_cml, p_cml)
  # FusioMR-s with uhp only
  niter = 20000
  K = length(b_exp)
  # starting values
  start_val = init_setup(niter, K, alpha_init = 1, beta_init = 0, sigma_gamma_init = 1, sigma_theta_init = 1)
  # estimation by R 
  # res = gibbs_single_uhp_only(niter, K, start_val$beta_tk, start_val$gamma_tk, start_val$theta_tk, start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk, b_out, b_exp, (se_out)^2, (se_exp)^2, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
  # estimation by Rcpp
  res = gibbs_seso_uhp_only_cpp(niter, K, start_val$beta_tk, start_val$theta_tk, start_val$gamma_tk, start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk, b_out, b_exp, (se_out)^2, (se_exp)^2, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
  if(all(!is.na(res))) {
    ids = (niter/2 + 1):niter
    bhat = mean(res$beta_tk[ids], na.rm = T)
    se_bhat = sd(res$beta_tk[ids], na.rm = T)
    pval = 2*(1-pnorm(abs(bhat)/se_bhat))
  } else { 
    bhat = se_bhat = pval = NA 
  }
  # summarize
  out_fusio = c(bhat, se_bhat, pval)
  output = c(gene, out_other_methods, out_fusio)
  names(output) = c("gene", "niv", "b_ivw", "se_ivw", "p_ivw", "b_egger", "se_egger", "p_egger", "b_cml", "se_cml", "p_cml", "b_fusio", "se_fusio", "p_fusio")
  return(output)
}

# function for sc only analysis
process_gene_sc_only <- function(dat1, gene_list1, j, pcut = 0.01, clump_kb = 10, clump_r2 = 0.1, b_gamma_prior_mean, b_theta_prior_mean) {
message(paste0("Processing gene-", j, ": ", gene_list1[j]))
tmp = dat1 %>% dplyr::filter(gene_symbol == gene_list1[j]) %>% dplyr::filter(p <= pcut)
if (nrow(tmp) > 0) {
# set hyper parameters S2
b_gamma_prior_mean = max(1e-3, mean(tmp$beta^2)-mean(tmp$se^2))
b_theta_prior_mean = max(1e-5, mean(tmp$b_out^2)-mean(tmp$se_out^2))
# Clumping
clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
# Run MR analysis if conditions are met
if (all(class(clumped_df) == c("data.table", "data.frame")) &&
nrow(clumped_df) > 1 && nrow(clumped_df) < 1000) {
b_exp = clumped_df$beta; se_exp = clumped_df$se; b_out = clumped_df$b_out; se_out = clumped_df$se_out
# set hyper parameters S1
K = nrow(clumped_df)
# b_gamma_prior_mean = max(1e-3, mean(b_exp^2)-mean(se_exp^2))
# b_theta_prior_mean = max(1e-5, mean(b_out^2)-mean(se_out^2))
a_gamma_prior = a_theta_prior = max(2, K/4)
b_gamma_prior = b_gamma_prior_mean * (a_gamma_prior-1)
b_theta_prior = b_theta_prior_mean * (a_theta_prior-1)
res = run_mr(gene_list1[j], b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
res = c(res, 'b_gamma_prior_mean' = b_gamma_prior_mean, 'b_theta_prior_mean' = b_theta_prior_mean)
return(res)
} else return(c(gene_list1[j], -1, rep(NA, 14))) # 0 or 1 snp after clumping
} else return(c(gene_list1[j], -2, rep(NA, 14))) # no snps w/ pval < pcut
}

# function for sc + bk
process_gene_sc_bk <- function(dat2, gene_list2, j, pcut = 0.01, clump_kb = 10, clump_r2 = 0.1) {
  message(paste0("Processing gene-", j, ": ", gene_list2[j]))
  # IVs = p < pcut | (p_bk < pcut & same sign)
  # tmp = dat2 %>% dplyr::filter(gene_symbol == gene_list2[j]) %>% dplyr::filter((p < pcut | (p_bk < pcut & sign(b_bk) == sign(beta))) & se > 0 & se_out > 0)
  # IVs = p < pcut | (p_bk < pcut)
  tmp = dat2 %>% dplyr::filter(gene_symbol == gene_list2[j]) %>% dplyr::filter((p < pcut | p_bk < pcut) & se > 0 & se_out > 0)
  if (nrow(tmp) > 0) {
    # Clumping
    clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
    # Run MR analysis if conditions are met
    if (all(class(clumped_df) == c("data.table", "data.frame")) &&
        nrow(clumped_df) >= 5 && nrow(clumped_df) < 1000) {
      res = run_mr(gene_list2[j], clumped_df$beta, clumped_df$se, clumped_df$b_out, clumped_df$se_out)
      return(res)
    }} else return(NULL) 
}


# function for sc + bk combine (add F_avg)
process_gene_sc_bk_combine <- function(dat3, gene_list3, j, pcut_sc = 0.01, pcut_bk = 0.001, clump_kb = 10, clump_r2 = 0.1, ncut_bk, b_gamma_prior_mean, b_theta_prior_mean) {
message(paste0("Processing gene-", j, ": ", gene_list3[j]))
# IVs = p < pcut | (p_bk < pcut in at least one tissue)
tmp = dat3 %>% dplyr::filter(gene_symbol == gene_list3[j]) 
nn = apply(tmp %>% dplyr::select(contains("p_")), 1, function(t) sum(t < pcut_bk, na.rm = T))
tmp = cbind(tmp, nn)
tmp = tmp %>% dplyr::filter((p < pcut_sc & nn >= ncut_bk))
if (nrow(tmp) > 0) {
# set hyper parameters S2
b_gamma_prior_mean = max(1e-3, mean(tmp$beta^2)-mean(tmp$se^2))
b_theta_prior_mean = max(1e-5, mean(tmp$b_out^2)-mean(tmp$se_out^2))
# Clumping
clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
# Run MR analysis if conditions are met
if (all(class(clumped_df) == c("data.table", "data.frame")) &&
nrow(clumped_df) >= 2 && nrow(clumped_df) < 1000) {
b_exp = clumped_df$beta; se_exp = clumped_df$se; b_out = clumped_df$b_out; se_out = clumped_df$se_out
# set hyper parameters S1
K = nrow(clumped_df)
# b_gamma_prior_mean = max(1e-3, mean(b_exp^2)-mean(se_exp^2))
# b_theta_prior_mean = max(1e-5, mean(b_out^2)-mean(se_out^2))
a_gamma_prior = a_theta_prior = max(2, K/4)
b_gamma_prior = b_gamma_prior_mean * (a_gamma_prior-1)
b_theta_prior = b_theta_prior_mean * (a_theta_prior-1)
res = run_mr(gene_list3[j], b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
F_avg = median(b_exp^2/se_exp^2)
res = c(res, 'b_gamma_prior_mean' = b_gamma_prior_mean, 'b_theta_prior_mean' = b_theta_prior_mean, 'F_avg' = F_avg)
return(res)
} else return(c(gene_list3[j], -1, rep(NA, 15))) # 0 or 1 snp after clumping
} else return(c(gene_list3[j], -2, rep(NA, 15))) # no snps w/ pval < pcut
}


assess_IV_strength <- function(dat1, gene_list1, j, pcut = 0.01, clump_kb = 10, clump_r2 = 0.1, nn = 192) {
message(paste0("Processing gene-", j, ": ", gene_list1[j]))
tmp = dat1 %>% dplyr::filter(gene_symbol == gene_list1[j]) %>% dplyr::filter(p <= pcut)
if (nrow(tmp) > 0) {
# Clumping
clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
# Run MR analysis if conditions are met
if (all(class(clumped_df) == c("data.table", "data.frame")) &&
nrow(clumped_df) > 1 && nrow(clumped_df) < 1000) {
b_exp = clumped_df$beta; se_exp = clumped_df$se
K = nrow(clumped_df)
R2 = sum(b_exp^2/(b_exp^2 + nn*se_exp^2))
if (R2 >= 1) F_stat = Inf else F_stat = R2/(1-R2)*(nn-K-1)/K
F_avg = median(b_exp^2/se_exp^2)
return(c(gene_list1[j], K, F_stat, F_avg))
} else return(c(gene_list1[j], -1, rep(NA,2))) # 0 or 1 snp after clumping
} else return(c(gene_list1[j], -2, rep(NA,2))) # no snps w/ pval < pcut
}

# function to run FusioMR seso only
run_mr_seso = function(gene, b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior) {
# mr.obj = MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
# IVW_f = MendelianRandomization::mr_ivw(mr.obj, model = 'fixed')
# b_ivw_fixed = IVW_f$Estimate; se_ivw_fixed = IVW_f$StdError; p_ivw_fixed = IVW_f$Pvalue
# Egger = tryCatch({MendelianRandomization::mr_egger(mr.obj)}, error = function(e) {NA})
# if(!is.null(Egger) & class(Egger)=='Egger') { 
# b_egger = Egger$Estimate; se_egger = Egger$StdError.Est; p_egger = Egger$Pvalue.Est
# } else b_egger = se_egger = p_egger = NA
# cml = tryCatch({MendelianRandomization::mr_cML(mr.obj, MA = TRUE, DP = FALSE, n = 100)}, error = function(e) {NA})
# if(!is.null(cml)) {  
# b_cml = cml$Estimate; se_cml = cml$StdError; p_cml = cml$Pvalue
# } else b_cml = se_cml = p_cml = NA
# out_other_methods = c(length(b_exp), b_ivw_fixed, se_ivw_fixed, p_ivw_fixed, b_egger, se_egger, p_egger, b_cml, se_cml, p_cml)
# FusioMR-s with uhp only
niter = 20000
K = length(b_exp)
# hyper parameters
# a_gamma_prior = b_gamma_prior = a_theta_prior = b_theta_prior = 0
# a_gamma_prior = a_theta_prior = 0
# b_gamma_prior = min(c(se_exp^2, se_out^2)/1)*(K/2) 
# b_theta_prior = min(c(se_exp^2, se_out^2)/1)*(K/2) 
# starting values
start_val = init_setup(niter, K, alpha_init = 1, beta_init = 0, sigma_gamma_init = 0.1, sigma_theta_init = 0.1)
# # estimation by R 
# res = gibbs_single_uhp_only(niter, K, start_val$beta_tk, start_val$gamma_tk, start_val$theta_tk, start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk, b_out, b_exp, (se_out)^2, (se_exp)^2, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
# estimation by Rcpp
res = gibbs_seso_uhp_only_cpp(niter, K, start_val$beta_tk, start_val$theta_tk, start_val$gamma_tk, start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk, b_out, b_exp, (se_out)^2, (se_exp)^2, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
if(all(!is.na(res))) {
ids = (niter/2 + 1):niter
bhat = mean(res$beta_tk[ids], na.rm = T)
se_bhat = sd(res$beta_tk[ids], na.rm = T)
pval = 2*(1-pnorm(abs(bhat)/se_bhat))
} else { 
bhat = se_bhat = pval = NA 
}
# summarize
out_fusio = c(bhat, se_bhat, pval)
output = c(gene, K, out_fusio)
names(output) = c("gene", "niv", "b_fusio", "se_fusio", "p_fusio")
return(output)
}

# function for bk only analysis
process_gene_bk_only <- function(dat1, gene_list1, j, pcut = 0.01, clump_kb = 10, clump_r2 = 0.1, b_gamma_prior_mean, b_theta_prior_mean) {
message(paste0("Processing gene-", j, ": ", gene_list1[j]))
tmp = dat1 %>% dplyr::filter(gene_name == gene_list1[j]) %>% dplyr::filter(p <= pcut)
if (nrow(tmp) > 0) {
# set hyper parameters S2
b_gamma_prior_mean = max(1e-3, mean(tmp$beta^2)-mean(tmp$se^2))
b_theta_prior_mean = max(1e-5, mean(tmp$b_out^2)-mean(tmp$se_out^2))
# Clumping
clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
# Run MR analysis if conditions are met
if (all(class(clumped_df) == c("data.table", "data.frame")) &&
nrow(clumped_df) > 1 && nrow(clumped_df) < 1000) {
b_exp = clumped_df$beta; se_exp = clumped_df$se; b_out = clumped_df$b_out; se_out = clumped_df$se_out
# set hyper parameters S1
K = nrow(clumped_df)
# b_gamma_prior_mean = max(1e-3, mean(b_exp^2)-mean(se_exp^2))
# b_theta_prior_mean = max(1e-5, mean(b_out^2)-mean(se_out^2))
a_gamma_prior = a_theta_prior = max(2, K/4)
b_gamma_prior = b_gamma_prior_mean * (a_gamma_prior-1)
b_theta_prior = b_theta_prior_mean * (a_theta_prior-1)
res = run_mr_seso(gene_list1[j], b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
res = c(res, 'b_gamma_prior_mean' = b_gamma_prior_mean, 'b_theta_prior_mean' = b_theta_prior_mean)
return(res)
} else return(c(gene_list1[j], -1, rep(NA, 5))) # 0 or 1 snp after clumping
} else return(c(gene_list1[j], -2, rep(NA, 5))) # no snps w/ pval < pcut
}

# function for sc only fujita analysis
process_gene_sc_only_fujita <- function(dat1, gene_list1, j, pcut = 0.01, clump_kb = 10, clump_r2 = 0.1, b_gamma_prior_mean, b_theta_prior_mean) {
message(paste0("Processing gene-", j, ": ", gene_list1[j]))
tmp = dat1 %>% dplyr::filter(gene_id == gene_list1[j]) %>% dplyr::filter(p <= pcut)
if (nrow(tmp) > 0) {
# set hyper parameters S2
b_gamma_prior_mean = max(1e-3, mean(tmp$beta^2)-mean(tmp$se^2))
b_theta_prior_mean = max(1e-5, mean(tmp$b_out^2)-mean(tmp$se_out^2))
# Clumping
clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
# Run MR analysis if conditions are met
if (all(class(clumped_df) == c("data.table", "data.frame")) &&
nrow(clumped_df) > 1 && nrow(clumped_df) < 1000) {
b_exp = clumped_df$beta; se_exp = clumped_df$se; b_out = clumped_df$b_out; se_out = clumped_df$se_out
# set hyper parameters S1
K = nrow(clumped_df)
# b_gamma_prior_mean = max(1e-3, mean(b_exp^2)-mean(se_exp^2))
# b_theta_prior_mean = max(1e-5, mean(b_out^2)-mean(se_out^2))
a_gamma_prior = a_theta_prior = max(2, K/4)
b_gamma_prior = b_gamma_prior_mean * (a_gamma_prior-1)
b_theta_prior = b_theta_prior_mean * (a_theta_prior-1)
F_avg = median(b_exp^2/se_exp^2)
res = run_mr_seso(gene_list1[j], b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
res = c(res, 'b_gamma_prior_mean' = b_gamma_prior_mean, 'b_theta_prior_mean' = b_theta_prior_mean, 'F_avg' = F_avg)
return(res)
} else return(c(gene_list1[j], -1, rep(NA, 6))) # 0 or 1 snp after clumping
} else return(c(gene_list1[j], -2, rep(NA, 6))) # no snps w/ pval < pcut
}




