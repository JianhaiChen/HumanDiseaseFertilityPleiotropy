## Figure 2 data prep — local significant-cells version, canonical antsyn (Ant=++/--).
## Displayed genes = IV-filtered keep_after_filter male dual genes per disease.
suppressMessages(library(data.table)); setwd("/Volumes/X10Pro/data/synergistic/mr_res")
## 1. displayed (gene,disease) = male dual, IV-kept; class from combo
kf<-fread("bonf050_broadmixed_recomb_block_shared_iv_filter_combos_with_annotation.csv")[keep_after_filter==TRUE & fitness=="male"]
kf[,cls:=fifelse(combo%in%c("++","--"),"ant","syn")]
disp<-kf[,.(disease,gene,conflict_type=combo,conflict_class=cls)]
cat("male dual (gene,disease) pairs:",nrow(disp)," genes:",uniqueN(disp$gene)," diseases:",uniqueN(disp$disease),"\n")
## 2. per-tissue z from sig_cells_long (male). use tissue_bonf standard (per-cell sig)
sc<-fread("dual_sig_standard_compare_sig_cells_long.csv")[sex=="male" & standard=="tissue_bonf",
   .(disease,gene,tissue,z_disease,z_fitness)]
## keep only displayed genes per disease
sc<-merge(sc,disp[,.(disease,gene)],by=c("disease","gene"))
cat("per-tissue cells for displayed genes:",nrow(sc)," (diseases w/ cells:",uniqueN(sc$disease),")\n")
## 3. per-disease summary + category + display name
info<-fread("disease_gwas_source_info.csv",select=c("disease","Disease","Category"))
cleanname<-function(s)gsub("[ /]+","_",trimws(gsub("\\s*\\([^)]*\\)","",s)))
info[,disp:=cleanname(Disease)]
mapcat<-function(c)fifelse(c=="neuro_psych","Neuropsychiatric",fifelse(c=="immune","Immune",
  fifelse(c%in%c("cardiovasc","metabolic"),"Cardiometabolic",fifelse(c=="cancer","Cancer",
  fifelse(c=="reproductive","Reproductive","Other")))))
info[,cat6:=mapcat(Category)]
dsum<-disp[,.(n_gene=uniqueN(gene),n_ant=uniqueN(gene[conflict_class=="ant"])),by=disease]
dsum[,ant_ratio:=n_ant/n_gene]
dsum<-merge(dsum,info[,.(disease,disp,cat6)],by="disease",all.x=TRUE)
dsum[is.na(disp),disp:=disease][is.na(cat6),cat6:="Other"]
## only diseases that actually have per-tissue cells & >=3 genes shown
have<-sc[,.(ncell=.N,ng=uniqueN(gene)),by=disease]
dsum<-merge(dsum,have,by="disease")[ng>=3]
cat("\ndiseases with >=3 displayed genes & cells:",nrow(dsum),"\n")
print(dsum[order(cat6,-ant_ratio),.(disease,disp,cat6,n_ant,n_gene,ant_ratio=round(ant_ratio,2),ng)])
saveRDS(list(disp=disp,sc=sc,dsum=dsum),"/Users/cjh/.claude/jobs/bb3be850/tmp/fig2_data.rds")
cat("\nsaved fig2_data.rds\n")
