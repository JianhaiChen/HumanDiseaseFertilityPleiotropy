suppressMessages(library(data.table))
MR<-"/Volumes/X10Pro/data/synergistic/mr_res"; AG<-"/Volumes/X10Pro/data/synergistic/mr_res/fusios_independent_iv/agg_traits"; strip<-function(x)sub("[.].*","",x); nk<-function(x)gsub("[^a-z0-9]","",tolower(x))
co<-fread(file.path(MR,"gene_coords_v19.tsv"),col.names=c("g","chr","pos"))[,g:=strip(g)]
isMHC<-function(chr,pos)chr==6 & pos>=25e6 & pos<=34e6
# IV过滤基因集(keep_after_filter)
kf<-unique(fread(file.path(MR,"bonf050_broadmixed_recomb_block_shared_iv_filter_combos_with_annotation.csv"))[keep_after_filter==TRUE,strip(gene)])
# --- FusioMR 权威: z=top_beta/top_se, MHC排除, |z|<9 ---
b1<-fread(file.path(MR,"gene_level_b1_disease_acat_summary.dropbox.csv"),select=c("gene","disease","top_beta","top_se"))[,g:=strip(gene)][is.finite(top_beta)&top_se>0][,z1:=top_beta/top_se]
b2<-fread(file.path(MR,"gene_level_b2_fitness_acat_summary.dropbox.csv"),select=c("gene","fitness","top_beta","top_se"))[,g:=strip(gene)][is.finite(top_beta)&top_se>0][,z2:=top_beta/top_se]
b1<-merge(b1,co,by="g")[!isMHC(chr,pos)&abs(z1)<9];b2<-merge(b2,co,by="g")[!isMHC(chr,pos)&abs(z2)<9]
mrRho<-function(genes){rbindlist(lapply(c("male","female"),function(sx){sxl<-fifelse(sx=="male","M","F");f<-b2[fitness==sx & (is.null(genes)|g%in%genes)]
  rbindlist(lapply(unique(b1$disease),function(dz){m<-merge(b1[disease==dz&(is.null(genes)|g%in%genes),.(g,z1)],f[,.(g,z2)],by="g");if(nrow(m)<30)return(NULL);data.table(dk=nk(dz),sex=sxl,rge=cor(m$z1,m$z2))}))}))}
mr_all<-dcast(mrRho(NULL),dk~sex,value.var="rge")
mr_iv<-dcast(mrRho(kf),dk~sex,value.var="rge")
# --- FusioS 权威方法: z=b_fusio/se_fusio, MHC排除, |z|<9 ---
rd<-function(t){f<-file.path(AG,paste0(t,".txt"));if(!file.exists(f))return(NULL);x<-fread(f);x[,g:=strip(gene)];x<-x[is.finite(b_fusio)&is.finite(se_fusio)&se_fusio>0][,z:=b_fusio/se_fusio];merge(unique(x[,.(g,z)],by="g"),co,by="g")[!isMHC(chr,pos)&abs(z)<9]}
fa<-rd("father");mo<-rd("mother");dzs<-setdiff(sub(".txt","",list.files(AG,"\\.txt$")),c("father","mother"))
fsRho<-function(genes){rbindlist(lapply(dzs,function(dz){d<-rd(dz);if(is.null(d))return(NULL)
  rf<-merge(d[is.null(genes)|g%in%genes,.(g,z1=z)],fa[,.(g,z2=z)],by="g");rm<-merge(d[is.null(genes)|g%in%genes,.(g,z1=z)],mo[,.(g,z2=z)],by="g")
  data.table(dk=nk(dz),M=if(nrow(rf)>=30)cor(rf$z1,rf$z2)else NA,F=if(nrow(rm)>=30)cor(rm$z1,rm$z2)else NA)}))}
fs_all<-fsRho(NULL);fs_iv<-fsRho(kf)
cat("=== 权威方法(z-score, MHC排除, |z|<9)ρ_GE 正比例 ===\n")
sm<-function(t,lab)cat(sprintf("  %-24s male 正%.0f%% | female 正%.0f%%\n",lab,100*mean(t$M>0,na.rm=T),100*mean(t$F>0,na.rm=T)))
sm(mr_all,"FusioMR 全基因");sm(mr_iv,"FusioMR IV过滤")
sm(fs_all,"FusioS 全基因");sm(fs_iv,"FusioS IV过滤")
# 对比 FusioMR全基因 vs FusioS全基因
c1<-merge(mr_all[,.(dk,mrM=M)],fs_all[,.(dk,fsM=M)],by="dk");c2<-merge(mr_all[,.(dk,mrF=F)],fs_all[,.(dk,fsF=F)],by="dk")
cat(sprintf("\nFusioMR vs FusioS (全基因, 权威法): male r=%.2f 符号一致%.0f%% | female r=%.2f 符号一致%.0f%%\n",
  cor(c1$mrM,c1$fsM,use="complete"),100*mean(sign(c1$mrM)==sign(c1$fsM),na.rm=T),cor(c2$mrF,c2$fsF,use="complete"),100*mean(sign(c2$mrF)==sign(c2$fsF),na.rm=T)))
fwrite(merge(mr_all,fs_all,by="dk",suffixes=c("_MR","_FS")),"/Users/cjh/Documents/GitHub/HumanDisease-FertilityPleiotropy/data/rhoge_auth_compare.csv")
