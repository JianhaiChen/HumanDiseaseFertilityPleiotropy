library(data.table)
library(tidyverse)
library(pbapply)
library(MendelianRandomization)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
source("code/functions.R")
#source("/gpfs/data/linchen-lab/Bowei/ad_educ/code/init_setup.R")
#source("/gpfs/data/linchen-lab/Bowei/ad_educ/code/label_flip.R")
#Rcpp::sourceCpp("/gpfs/data/linchen-lab/Bowei/ad_educ/code/gibbs_joint_rcpp_nopr.cpp", cacheDir ="/gpfs/data/linchen-lab/Bowei/")
#Rcpp::sourceCpp("/gpfs/data/linchen-lab/Bowei/ad_educ/code/gibbs_rcpp_pp.cpp", cacheDir ="/gpfs/data/linchen-lab/Bowei/")
source("/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code/init_setup_seso.R")
#Rcpp::sourceCpp("/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code/gibbs_seso_uhp_only.cpp", cacheDir ="/gpfs/data/linchen-lab/Bowei/")
source("/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code/gibbs_single_uhp_only.R")

args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

pcut = 0.01
clump_kb = 50; clump_r2 = 0.01
#clump_param = list(clump_kb = c(10,20), clump_r2 = c(0.1,0.1))
ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
diseases_short = "ad"

#for (k in 1:2) {
#clump_kb = clump_param$clump_kb[k]
#clump_r2 = clump_param$clump_r2[k]
#print(paste0("clump param ", k, ", kb=", clump_kb, ", r2=", clump_r2))
for (i in 1:length(ctypes_short)) {
print(paste0("ctype-", i, ", ", ctypes_short[i]))
dat = fread(paste0('harmo_out/eqtl_', diseases_short, '_', ctypes_short[i], "_", chr, '.txt.gz')) 
sumtbl1 = NULL

# IVs are selected from sc only
dat1 = dat %>% dplyr::select(gene_symbol:se, b_out, se_out) %>% dplyr::filter(se != Inf, se_out != Inf, se > 0, se_out > 0) %>% drop_na() 
# aa = dat1 %>% group_by(gene_name) %>% summarise(n = n()); summary(aa$n)
gene_list1 = unique(dat1$gene_symbol)
# set hyper parameters S3
b_gamma_prior_mean = max(1e-3, mean(dat1$beta^2)-mean(dat1$se^2))
b_theta_prior_mean = max(1e-5, mean(dat1$b_out^2)-mean(dat1$se_out^2))
sumtbl1 = do.call(rbind, pblapply(seq_along(gene_list1), function(t) process_gene_sc_only(dat1, gene_list1, t, pcut, clump_kb, clump_r2, b_gamma_prior_mean, b_theta_prior_mean)))
fwrite(sumtbl1, paste0("sumtbl_sconly_uhp/sconlys2_", diseases_short, "_", pcut, "_", clump_kb, "_", clump_r2, "_", ctypes_short[i], "_", chr, '.txt'))
} # end loop cell type
#} # end loop clumping parameters


