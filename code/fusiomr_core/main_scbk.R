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
tissue = as.character(args[1])
chr = as.numeric(args[2])
print(tissue)
print(chr)

brain_tissues = c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra")

pcut = 0.001
# clump_kb = 10; clump_r2 = 0.1
clump_param = list(clump_kb = c(50,50), clump_r2 = c(0.1,0.01))
# ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
ctypes_short = c("end")
diseases_short = "ad"

for (k in 1:2) {
clump_kb = clump_param$clump_kb[k]
clump_r2 = clump_param$clump_r2[k]
print(paste0("clump param ", k, ", kb=", clump_kb, ", r2=", clump_r2))
for (i in 1:length(ctypes_short)) {
print(paste0("ctype-", i, ", ", ctypes_short[i]))
dat = fread(paste0('harmo_out/eqtl_', diseases_short, '_', ctypes_short[i], "_", chr, '.txt.gz'))
sumtbl2 = NULL
dat2 = dat %>% dplyr::select(gene_symbol:se, contains(tissue), b_out, se_out) %>% drop_na()
colnames(dat2)[9:11] = paste(c('p','b','se'), 'bk', sep = '_')
# aa = dat2 %>% group_by(gene_name) %>% summarise(n = n()); summary(aa$n)
gene_list2 = unique(dat2$gene_symbol)
sumtbl2 = do.call(rbind, pblapply(seq_along(gene_list2), function(t) process_gene_sc_bk(dat2, gene_list2, t, pcut, clump_kb, clump_r2)))
fwrite(sumtbl2, paste0("sumtbl_uhp/scbk_", diseases_short, "_", pcut, "_", clump_kb, "_", clump_r2, "_", ctypes_short[i], "_", chr, "_", tissue, ".txt"))
} # end loop cell type
} # end loop clumping parameters


