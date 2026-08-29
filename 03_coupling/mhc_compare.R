suppressMessages(library(data.table))
MR<-"${MR}"; strip<-function(x)sub("[.].*","",x); nk<-function(x)gsub("[^a-z0-9]","",tolower(x))
co<-fread(file.path(MR,"gene_coords_v19.tsv"),col.names=c("g","chr","pos"))[,g:=strip(g)]
isM<-function(chr,pos)chr==6&pos>=25e6&pos<=34e6
b1<-fread(file.path(MR,"gene_level_b1_disease_acat_summary.dropbox.csv"),select=c("gene","disease","top_beta","top_se"))[,g:=strip(gene)][is.finite(top_beta)&top_se>0][,z1:=top_beta/top_se];b1<-merge(b1,co,by="g")[abs(z1)<9]
b2<-fread(file.path(MR,"gene_level_b2_fitness_acat_summary.dropbox.csv"),select=c("gene","fitness","top_beta","top_se"))[,g:=strip(gene)][is.finite(top_beta)&top_se>0][,z2:=top_beta/top_se];b2<-merge(b2,co,by="g")[abs(z2)<9]
b1[,mhc:=isM(chr,pos)]; b2[,mhc:=isM(chr,pos)]
rge<-function(useMHC){rbindlist(lapply(unique(b1$disease),function(dz){rbindlist(lapply(c("male","female"),function(sx){
  x1<-b1[disease==dz]; x2<-b2[fitness==sx]; if(!useMHC){x1<-x1[!mhc==T];x2<-x2[mhc==F]}
  m<-merge(x1[,.(g,z1)],x2[,.(g,z2)],by="g");if(nrow(m)<30)return(NULL)
  data.table(dk=nk(dz),sex=sx,rge=cor(m$z1,m$z2),ng=nrow(m))}))}))}
noM<-rge(FALSE)[,.(dk,sex,rge_noMHC=rge)]; wM<-rge(TRUE)[,.(dk,sex,rge_MHC=rge,ng=ng)]
d<-merge(noM,wM,by=c("dk","sex"))
# 免疫病标记
imm<-c("graves","t1d","ulcerativecolitis","crohns","hayfever","primarybiliary","rheumatoid","sle","psoriasis","multiplesclerosis","celiac","atopiderm","eczema","asthma","alopecia","vitiligo","ankylosing","lupus","hashimoto","hypothyroid","hodgkin","lymphoma")
d[,immune:=grepl(paste(imm,collapse="|"),dk)]
d[,delta:=rge_MHC-rge_noMHC]
cat("=== MHC 纳入的影响: 免疫病 vs 非免疫病 |Δρ_GE| ===\n")
cat(sprintf("  免疫病 中位|Δ|=%.4f (n=%d)  非免疫病 中位|Δ|=%.4f (n=%d)  Wilcox P=%.2g\n",
  median(abs(d[immune==T]$delta)),nrow(d[immune==T]),median(abs(d[immune==F]$delta)),nrow(d[immune==F]),
  wilcox.test(abs(d[immune==T]$delta),abs(d[immune==F]$delta))$p.value))
cat(sprintf("  免疫病 中位Δ(带符号)=%+.4f  非免疫病=%+.4f  (>0=纳入后升)\n",median(d[immune==T]$delta),median(d[immune==F]$delta)))
cat("\n=== 免疫病 MHC 前后 (男女) ===\n")
print(d[immune==T,.(dk,sex,noMHC=round(rge_noMHC,4),withMHC=round(rge_MHC,4),delta=round(delta,4))][order(-abs(delta))][1:16])
cat("\n=== 整体相关是否稳定(非免疫病) ===\n")
ni<-d[immune==F];cat(sprintf("  cor(noMHC, withMHC) 非免疫病 = %.4f\n",cor(ni$rge_noMHC,ni$rge_MHC)))
cat(sprintf("  全部病 cor = %.4f\n",cor(d$rge_noMHC,d$rge_MHC)))
fwrite(d,"mhc_compare_table.csv")
