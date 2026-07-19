library(data.table); library(tidyverse)
args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
diseases = c("dcfdx_ad.assoc.logistic")
diseases_short = c("ad")

for (k in 1:length(diseases)) {
print(paste0(k, ", ", diseases[k]))
out = fread(paste0("../ad_educ/ad_data/", diseases[k], ".gz_all.txt")) %>% dplyr::select(-pos) %>% dplyr::rename(b_out=beta, se_out=se)

for (i in 1:length(ctypes_short)) {
print(paste0(i, ", ", ctypes_short[i]))
tmp = fread(paste0("harmo/eqtl_", ctypes_short[i], "_", chr, '.txt.gz')) %>% left_join(out, by = c('rsid'='rsid')) %>% mutate(align = ifelse(effect_allele==EFF & other_allele==REF, 1, ifelse(effect_allele==REF & other_allele==EFF, -1, NA))) %>% dplyr::filter(!(is.na(align) & !(is.na(b_out)))) %>% mutate(b_out = b_out*align) %>% drop_na(b_out, se_out) %>% dplyr::select(-REF, -EFF, -align, -chr)
fwrite(tmp, paste0('harmo_out/eqtl_', diseases_short[k], '_', ctypes_short[i], "_", chr, '.txt'))
cmd_zip = paste0('gzip harmo_out/eqtl_', diseases_short[k], '_', ctypes_short[i], "_", chr, '.txt')
system(cmd_zip)
} # end loop cell type
} # end loop disease

