library(data.table); library(tidyverse)

args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

snps = fread("bryois/snp_pos.txt.gz", select = c("SNP", "effect_allele", "other_allele"))
brain_tissues = c("Brain_Amygdala",
"Brain_Anterior_cingulate_cortex_BA24",
"Brain_Caudate_basal_ganglia",
"Brain_Cerebellar_Hemisphere",
"Brain_Cerebellum",
"Brain_Cortex",
"Brain_Frontal_Cortex_BA9",
"Brain_Hippocampus",
"Brain_Hypothalamus",
"Brain_Nucleus_accumbens_basal_ganglia",
"Brain_Putamen_basal_ganglia",
"Brain_Spinal_cord_cervical_c-1",
"Brain_Substantia_nigra")

#ctypes = c("Endothelial.cells", "Astrocytes", "Excitatory.neurons", "Inhibitory.neurons", "Microglia", "OPCs...COPs", "Oligodendrocytes", "Pericytes")
#ctypes_short = c("end", "ast", "exc", "inh", "mic", "opc", "oli", "per")
ctypes = c("Excitatory.neurons", "Inhibitory.neurons", "Microglia", "OPCs...COPs", "Oligodendrocytes", "Pericytes")
ctypes_short = c("exc", "inh", "mic", "opc", "oli", "per")

for (i in 1:length(ctypes)) {
print(paste0(i, ", ", ctypes[i]))
# colnames: Gene_id, SNP_id, Distance to TSS, Nominal p-value, Beta
tmp = fread(paste0("../single_cell_eQTL/bryois/", ctypes[i], ".", chr, ".gz")) %>% inner_join(snps, by = c("V2"="SNP")) %>% drop_na() %>% separate(V1, into = c("gene_symbol", "gene_name"), sep = "_") %>% dplyr::rename(rsid=V2, p=V4, beta=V5) %>% mutate(se = abs(beta/qnorm(p/2, lower.tail = F))) %>% dplyr::select(-V3)
for (j in 1:length(brain_tissues)) {
print(paste0(j, ", ", brain_tissues[j]))
eqtl_bk = fread(paste0("gtex/", brain_tissues[j], "_", chr, '.txt.gz')) %>% separate(gene_id, into = c("gene_name", "transcript"), sep = "\\.") 
# need to handle 1-to-many map
tmp = tmp %>% left_join(eqtl_bk, by = c("gene_name", "rsid")) %>% mutate(align = ifelse(effect_allele==alt & other_allele==ref, 1, ifelse(effect_allele==ref & other_allele==alt, -1, NA))) %>% dplyr::filter(!(is.na(align) & !(is.na(slope)))) %>% mutate(slope = slope*align) %>% dplyr::select(-ref, -alt, -transcript, -align, -chr)
colnames(tmp)[(ncol(tmp)-2):ncol(tmp)] = paste(c('p','b','se'), brain_tissues[j], sep = '_')
} # end loop tissues
fwrite(tmp, paste0('harmo/eqtl_', ctypes_short[i], '_', chr, '.txt'))
cmd_zip = paste0("gzip harmo/eqtl_", ctypes_short[i], '_', chr, '.txt')
system(cmd_zip)
} # end loop cell types




