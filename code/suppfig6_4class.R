# FigS_rg_4class_AB (recent mode): rows A=r_g(LDSC) / B=rho_GE(MR); columns = 4 effect-direction classes.
# Each point = one disease-sex; y = % of that disease's pleiotropic genes in the class (>=10 pleiotropic genes).
# Sex: male solid/triangle, female dashed/cross; weighted LM (weight = #pleiotropic genes); r,P per sex per panel.
suppressMessages({library(data.table);library(ggplot2);library(patchwork)})
setwd("/Volumes/X10Pro/data/synergistic/mr_res")
CFPATH<-"/Volumes/X10Pro/data/synergistic/mr_res/session_2026-07-09_iv_filter/colleague_shared_iv_filter.csv"
BS<-8;th<-theme_classic(base_size=BS)+theme(axis.title=element_text(size=BS-1),axis.text=element_text(size=BS-2,color="grey20"),legend.text=element_text(size=BS-1),panel.grid=element_blank())
co<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[chr%in%as.character(1:22)][,.(gene,chr=as.integer(chr),pos=as.numeric(pos))]
isMHC<-function(chr,pos)chr==6&pos>25e6&pos<34e6
b1<-fread("gene_level_b1_disease_acat_summary.dropbox.csv",select=c("gene","disease","top_beta","top_se"))
b2<-fread("gene_level_b2_fitness_acat_summary.dropbox.csv",select=c("gene","fitness","top_beta","top_se"))[!is.na(fitness)]
## ---- disease-fertility correlation per disease-sex: rho-GE (all genes) + r_g (LDSC) ----
za<-merge(b1[is.finite(top_beta)&top_se>0,.(gene,disease,z1=top_beta/top_se)],co,by="gene")[!isMHC(chr,pos)&abs(z1)<9]
zb<-merge(b2[is.finite(top_beta)&top_se>0,.(gene,sex=fifelse(fitness=="male","M","F"),z2=top_beta/top_se)],co,by="gene")[!isMHC(chr,pos)&abs(z2)<9]
rge<-merge(za,zb[,.(gene,sex,z2)],by="gene",allow.cartesian=T)[,.(corr=cor(z1,z2),ng=.N),by=.(disease,sex)][ng>=30][,method:="expression effect (rho-GE, MR)"]
Lr<-readLines("dz_rg_97_ALLESTIMABLE.txt");vv<-function(p,t){m<-grep(paste0("^",t,"="),p,value=T);if(!length(m))return(NA);suppressWarnings(as.numeric(sub(paste0("^",t,"=[[:space:]]*(-?[0-9.eE-]+).*"),"\\1",m)))}
rgt<-rbindlist(lapply(Lr,function(x){p<-strsplit(x,"\t")[[1]];data.table(disease=sub("_$","",p[1]),M=vv(p,"father"),F=vv(p,"mother"))}))[!is.na(disease)&disease!=""&is.finite(M)&is.finite(F)]
rgl<-melt(rgt,id.vars="disease",variable.name="sex",value.name="corr")[,`:=`(sex=as.character(sex),method="genome-wide SNPs (r_g, LDSC)")]
CO<-rbind(rge[,.(disease,sex,corr,method)],rgl[,.(disease,sex,corr,method)])
## ---- four-class % per disease-sex, from IV-filtered pleiotropic set; keep >=10 pleiotropic genes ----
CF<-fread(CFPATH)[keep_after_filter==TRUE];CF[,sex:=fifelse(fitness=="male","M","F")]
rate<-CF[,.N,by=.(disease,sex,combo)][,tot:=sum(N),by=.(disease,sex)][tot>=10][,pct:=100*N/tot]
rate<-melt(dcast(rate,disease+sex+tot~combo,value.var="pct",fill=0),id.vars=c("disease","sex","tot"),variable.name="combo",value.name="pct")
db<-merge(rate,CO,by=c("disease","sex"),allow.cartesian=T)
db[,method:=factor(method,levels=c("genome-wide SNPs (r_g, LDSC)","expression effect (rho-GE, MR)"))]
db[,combo:=factor(combo,levels=c("++","-+","--","+-"),labels=c("d+ f+","d- f+","d- f-","d+ f-"))]
db[,sexL:=factor(fifelse(sex=="M","Male","Female"),levels=c("Male","Female"))]
## weighted r & P per method x class x sex (weight = # pleiotropic genes)
rC<-db[,{w<-tot;rw<-tryCatch(cov.wt(cbind(corr,pct),wt=w,cor=TRUE)$cor[1,2],error=function(e)NA_real_);pw<-tryCatch(coef(summary(lm(pct~corr,weights=w)))[2,4],error=function(e)NA_real_);.(r=rw,p=pw)},by=.(method,combo,sexL)]
rC[,lab:=sprintf("%s r=%.2f, P=%.0e",ifelse(sexL=="Male","M","F"),r,p)]
print(rC[order(method,combo,sexL)])
cols<-c(Male="#2C5985",Female="#C0392B")
mkrow<-function(meth,xl){d<-db[method==meth];rc<-rC[method==meth]
  ggplot(d,aes(corr,pct,color=sexL,fill=sexL,linetype=sexL,shape=sexL))+
    geom_point(size=.5,alpha=.45)+geom_smooth(aes(weight=tot),method="lm",formula=y~x,se=T,linewidth=.6,alpha=.12)+
    facet_wrap(~combo,nrow=1)+
    geom_text(data=rc[sexL=="Male"],aes(x=-Inf,y=Inf,label=lab),color=cols[["Male"]],hjust=-.06,vjust=1.4,size=1.8,inherit.aes=FALSE)+
    geom_text(data=rc[sexL=="Female"],aes(x=-Inf,y=Inf,label=lab),color=cols[["Female"]],hjust=-.06,vjust=2.9,size=1.8,inherit.aes=FALSE)+
    scale_color_manual(values=cols,name=NULL)+scale_fill_manual(values=cols,name=NULL)+scale_linetype_manual(values=c(Male="solid",Female="22"),name=NULL)+scale_shape_manual(values=c(Male=17,Female=4),name=NULL)+
    labs(x=xl,y="% pleiotropic genes")+th+theme(strip.text=element_text(face="bold",size=BS-1),strip.background=element_rect(fill="grey92",color=NA))}
pR<-mkrow("genome-wide SNPs (r_g, LDSC)",expression(italic(r)[g]*" (LDSC, genome-wide SNPs)"))
pG<-mkrow("expression effect (rho-GE, MR)",expression(rho[GE]*" (MR, expression effect)"))
fig<-pR/pG+plot_layout(guides="collect")+plot_annotation(tag_levels="A")&theme(legend.position="bottom",plot.tag=element_text(face="bold",size=12))
OUT<-"/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait"
ggsave(file.path(OUT,"FigS_rg_4class_AB.pdf"),fig,width=180,height=118,units="mm",device=cairo_pdf)
ggsave(file.path(OUT,"FigS_rg_4class_AB.png"),fig,width=180,height=118,units="mm",dpi=190)
cat("saved; disease-sex strata (>=10 genes):",uniqueN(db[,paste(disease,sex)]),"\n")
