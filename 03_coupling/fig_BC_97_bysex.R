# Fig3 B+C UPDATED 2026-07-12: 97 diseases (all-estimable r_g incl FinnGen rescues), SEX-SPECIFIC (no M/F averaging).
suppressMessages({library(data.table);library(ggplot2);library(patchwork);library(readxl)})
setwd("${MR}")
# STABLE paths (no ephemeral scratchpad): recomb maps in mr_res/recomb; colleague filter in session dir
SC<-"${MR}"                                      # recomb maps: SC/recomb/
CFPATH<-"${MR}/session_2026-07-09_iv_filter/colleague_shared_iv_filter.csv"
MY<-"${PROJ}"
BS<-8; th<-theme_classic(base_size=BS)+theme(axis.title=element_text(size=BS),axis.text=element_text(size=BS-1,color="grey20"),legend.text=element_text(size=BS-1),legend.key.size=unit(3.2,"mm"),panel.grid=element_blank())
co<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[chr%in%as.character(1:22)][,.(gene,chr=as.integer(chr),pos=as.numeric(pos))]
isMHC<-function(chr,pos) chr==6 & pos>25e6 & pos<34e6
b1<-fread("gene_level_b1_disease_acat_summary.dropbox.csv",select=c("gene","disease","top_beta","top_se","p_acat"))
b2<-fread("gene_level_b2_fitness_acat_summary.dropbox.csv",select=c("gene","fitness","top_beta","top_se","p_acat"))[!is.na(fitness)]
## ---- rho-GE: ALL genes, MHC excluded, |z|<9 cap ----
za<-merge(b1[is.finite(top_beta)&top_se>0,.(gene,disease,z1=top_beta/top_se)],co,by="gene")[!isMHC(chr,pos)&abs(z1)<9]
zb<-merge(b2[is.finite(top_beta)&top_se>0,.(gene,sex=fifelse(fitness=="male","M","F"),z2=top_beta/top_se)],co,by="gene")[!isMHC(chr,pos)&abs(z2)<9]
rge<-merge(za[,.(gene,disease,z1)],zb[,.(gene,sex,z2)],by="gene",allow.cartesian=T)[,.(rge=cor(z1,z2),ng=.N),by=.(disease,sex)][ng>=30]
rho<-dcast(rge,disease~sex,value.var="rge")
## ---- r_g (LDSC): 97 all-estimable, sex-specific father(M)/mother(F) ----
Lr<-readLines("dz_rg_97_ALLESTIMABLE.txt");vv<-function(p,t){m<-grep(paste0("^",t,"="),p,value=T);if(!length(m))return(c(NA,NA));c(suppressWarnings(as.numeric(sub(paste0("^",t,"=[[:space:]]*(-?[0-9.eE-]+).*"),"\\1",m))),suppressWarnings(as.numeric(sub(".*\\(([0-9.eE-]+)\\).*","\\1",m))))}
rgt<-rbindlist(lapply(Lr,function(x){p<-strsplit(x,"\t")[[1]];fa<-vv(p,"father");mo<-vv(p,"mother");data.table(disease=sub("_$","",p[1]),M=fa[1],se_fa=fa[2],F=mo[1],se_mo=mo[2])}))[!is.na(disease)&disease!=""&is.finite(M)&is.finite(F)][,`:=`(zM=M/se_fa,zF=F/se_mo)]
## ---- antagonism rate ----
CF<-fread(CFPATH)[keep_after_filter==TRUE]
rate<-CF[,.(n=.N,pct=100*mean(combo%in%c("++","--"))),by=.(disease,sex=fifelse(fitness=="male","M","F"))][n>=1]  # all 97; low-n de-emphasized by size & n-weighted regression
## ---- category map ----
nk<-function(x)gsub("[^a-z0-9]","",tolower(x));st1<-as.data.table(read_excel("Supplementary_Table_1.xlsx",sheet=1,skip=2));setnames(st1,2:4,c("Trait","Type","Category"));st1<-st1[grepl("Disease",Type,ignore.case=T),.(k=nk(Trait),Category)]
brc<-fread("Supplementary_Table_1_corrected_v2.csv",select=c("b1_disease_id","Suppl_Trait"));brc[,k:=nk(Suppl_Trait)];brc<-merge(brc,st1,by="k",all.x=T)
map6<-function(c)fifelse(c=="Neuropsychiatric","Neuropsychiatric",fifelse(grepl("Immune",c),"Immune",fifelse(c%in%c("Cardiovascular","Metabolic"),"Cardiometabolic",fifelse(c=="Cancers","Cancer",fifelse(c=="Reproductive","Reprod","Other")))))
brc[,cat6:=map6(Category)];dcat<-unique(brc[,.(k2=nk(b1_disease_id),cat6)])
## ================= Panel A (was B): correlation vs %antagonistic, by sex =================
CO<-rbind(rge[,.(disease,sex,corr=rge,method="gene expression (rho-GE, MR)")],melt(rgt[,.(disease,M,F)],id.vars="disease",variable.name="sex",value.name="corr")[is.finite(corr)][,`:=`(sex=as.character(sex),method="genome-wide SNPs (r_g, LDSC)")])
db<-merge(rate,CO,by=c("disease","sex"),allow.cartesian=T);db[,method:=factor(method,levels=c("gene expression (rho-GE, MR)","genome-wide SNPs (r_g, LDSC)"))];db[,sex:=factor(sex,levels=c("M","F"),labels=c("Male","Female"))]
rB<-db[,{w<-n;rw<-suppressWarnings(cov.wt(cbind(corr,pct),wt=w,cor=TRUE)$cor[1,2]);pw<-tryCatch(coef(summary(lm(pct~corr,weights=w)))[2,4],error=function(e)NA);.(r=round(rw,2),p=pw)},by=.(method,sex)];rB[,lab:=sprintf("%s: r=%.2f, P=%.1e",ifelse(sex=="Male","M","F"),r,p)]
## A-left: rho_GE vs %antagonistic (by sex); A-right: rho_GE vs r_g convergence (replaces old redundant r_g-vs-%antag facet -> moved to FigS validation)
dbL<-db[method=="gene expression (rho-GE, MR)"];rBL<-rB[method=="gene expression (rho-GE, MR)"]
pA1<-ggplot(dbL,aes(corr,pct,color=sex,fill=sex))+geom_point(aes(size=n),alpha=.5)+geom_smooth(aes(weight=n,linetype=sex),method="lm",formula=y~x,se=T,linewidth=.8,alpha=.12)+
  geom_text(data=rBL[sex=="Male"],aes(x=-Inf,y=Inf,label=lab),color="#2C5985",hjust=-.06,vjust=1.5,size=2.3,inherit.aes=F)+geom_text(data=rBL[sex=="Female"],aes(x=-Inf,y=Inf,label=lab),color="#C0392B",hjust=-.06,vjust=3,size=2.3,inherit.aes=F)+
  scale_size_continuous(range=c(.35,2),name="indep. loci",breaks=c(1,5,20,50))+
  scale_color_manual(values=c(Male="#2C5985",Female="#C0392B"),name=NULL)+scale_fill_manual(values=c(Male="#2C5985",Female="#C0392B"),name=NULL)+scale_linetype_manual(values=c(Male="solid",Female="22"),name=NULL)+guides(fill="none")+
  labs(x=expression("disease-fertility "*rho[GE]*" (gene expression)"),y="% antagonistic (indep. loci)")+th+theme(legend.position=c(.995,.02),legend.justification=c(1,0),legend.background=element_rect(fill=alpha("white",.7),color=NA),legend.box="horizontal")
conv<-merge(rgt[,.(disease,rgM=M,rgF=F)],rho[,.(disease,geM=M,geF=F)],by="disease")
convL<-rbind(conv[,.(sex="Male",rg=rgM,ge=geM)],conv[,.(sex="Female",rg=rgF,ge=geF)])[is.finite(rg)&is.finite(ge)][,sex:=factor(sex,levels=c("Male","Female"))]
rC<-convL[,{ct<-suppressWarnings(cor.test(rg,ge));.(lab=sprintf("%s: r=%.2f, P=%.1e",sex,ct$estimate,ct$p.value))},by=sex]
pA2<-ggplot(convL,aes(rg,ge,color=sex))+geom_hline(yintercept=0,linetype=2,color="grey85")+geom_vline(xintercept=0,linetype=2,color="grey85")+
  geom_point(size=1,alpha=.6)+geom_smooth(aes(linetype=sex),method="lm",formula=y~x,se=F,linewidth=.7)+
  geom_text(data=rC[sex=="Male"],aes(x=-Inf,y=Inf,label=lab),color="#2C5985",hjust=-.08,vjust=1.5,size=2.3,inherit.aes=F)+geom_text(data=rC[sex=="Female"],aes(x=-Inf,y=Inf,label=lab),color="#C0392B",hjust=-.08,vjust=3,size=2.3,inherit.aes=F)+
  scale_color_manual(values=c(Male="#2C5985",Female="#C0392B"),name=NULL)+scale_linetype_manual(values=c(Male="solid",Female="22"),name=NULL)+
  labs(x=expression(italic(r)[g]*" (LDSC, genome-wide SNPs)"),y=expression(rho[GE]*" (MR, gene expression)"))+th+theme(legend.position="none")
pB<-pA1|pA2
## ================= Panel B (was C): per-disease r_g & rho-GE, SEX-FACETED =================
aC<-merge(rgt[,.(disease,rg_M=M,rg_F=F,zM,zF)],rho[,.(disease,ge_M=M,ge_F=F)],by="disease",all.x=T)
aC[,k:=nk(disease)];aC<-merge(aC,brc[,.(k2=nk(b1_disease_id),cat6,Suppl_Trait)],by.x="k",by.y="k2",all.x=T)
aC[is.na(cat6),cat6:="Other"];aC[is.na(Suppl_Trait),Suppl_Trait:=disease];aC[,Suppl_Trait:=gsub("_"," ",Suppl_Trait)]
aC[,ord:=rowMeans(cbind(rg_M,rg_F),na.rm=T)]
rgcat<-aC[,.(m=mean(ord,na.rm=T)),by=cat6][order(-m)];aC[,cat6:=factor(cat6,levels=c(setdiff(rgcat$cat6,"Other"),"Other"))]
aC<-aC[order(cat6,ord)];aC[,xpos:=1:.N];bd<-aC[,.(mid=mean(xpos),x1=max(xpos)),by=cat6][order(cat6)];seps<-head(bd$x1,-1)+.5
dC<-rbindlist(list(
  aC[,.(Suppl_Trait,cat6,xpos,sex="Male",  measure="r_g",   corr=rg_M,sig=abs(zM)>1.96)],
  aC[,.(Suppl_Trait,cat6,xpos,sex="Male",  measure="rho_GE",corr=ge_M,sig=abs(zM)>1.96)],
  aC[,.(Suppl_Trait,cat6,xpos,sex="Female",measure="r_g",   corr=rg_F,sig=abs(zF)>1.96)],
  aC[,.(Suppl_Trait,cat6,xpos,sex="Female",measure="rho_GE",corr=ge_F,sig=abs(zF)>1.96)]))
dC[,sex:=factor(sex,levels=c("Male","Female"))]
dC[,corr_s:=corr/max(abs(corr),na.rm=T),by=measure];dC[,xplot:=xpos+ifelse(measure=="r_g",-.22,.22)]
dC[,grp:=fifelse(!sig|is.na(sig),"n.s.",fifelse(measure=="r_g","genome-wide SNPs (r_g, LDSC)","gene expression (rho-GE)"))]
dC[,grp:=factor(grp,levels=c("genome-wide SNPs (r_g, LDSC)","gene expression (rho-GE)","n.s."))]
cols<-c("genome-wide SNPs (r_g, LDSC)"="#2C5985","gene expression (rho-GE)"="#C0392B","n.s."="grey65");shp<-c("genome-wide SNPs (r_g, LDSC)"=16,"gene expression (rho-GE)"=16,"n.s."=1)
pC<-ggplot(dC,aes(xplot,corr_s,color=grp,shape=grp))+geom_hline(yintercept=0,color="grey55",linewidth=.3)+geom_vline(xintercept=seps,linetype="dashed",color="grey70",linewidth=.3)+geom_linerange(aes(ymin=0,ymax=corr_s),linewidth=.3)+geom_point(size=.9,stroke=.4)+
  facet_grid(sex~.)+geom_text(data=rbind(cbind(bd,sex="Male"),cbind(bd,sex="Female"))[,sex:=factor(sex,levels=c("Male","Female"))],aes(x=mid,y=.92,label=cat6),inherit.aes=FALSE,fontface="bold",size=1.9)+scale_color_manual(values=cols,name=NULL)+scale_shape_manual(values=shp,name=NULL)+
  scale_x_continuous(breaks=aC$xpos,labels=aC$Suppl_Trait,expand=expansion(add=.6))+coord_cartesian(ylim=c(-1.03,1.03),clip="off")+labs(x=NULL,y="disease-fertility correlation (scaled)")+
  theme_bw(7)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=4),panel.grid=element_blank(),strip.text=element_text(face="bold",size=8),strip.background=element_rect(fill="grey92",color=NA),legend.position="top",legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=5.6),legend.margin=margin(0,0,4,0),plot.margin=margin(1,4,1,1))
fig<-pB/pC+plot_layout(heights=c(1,1.5))+plot_annotation(tag_levels=list(c("A","B","C")))&theme(plot.tag=element_text(face="bold",size=12))
ggsave(file.path(MY,"Fig3_BC_97_bysex.pdf"),fig,width=182,height=185,units="mm",device=cairo_pdf)
ggsave(file.path(MY,"Fig3_BC_97_bysex.png"),fig,width=182,height=185,units="mm",dpi=185)
## ALSO write canonical main-text Fig 3 name (replaces stale 87-disease averaged version)
ggsave(file.path(MY,"Fig_BC_final.pdf"),fig,width=182,height=185,units="mm",device=cairo_pdf)
ggsave(file.path(MY,"Fig_BC_final.png"),fig,width=182,height=185,units="mm",dpi=185)
ggsave("${MR}/Fig_BC_final.pdf",fig,width=182,height=185,units="mm",device=cairo_pdf)
## sex-specific concordance stats
for(s in c("M","F")){cc<-aC[is.finite(get(paste0("rg_",s)))&is.finite(get(paste0("ge_",s)))];r<-cor(cc[[paste0("rg_",s)]],cc[[paste0("ge_",s)]]);cat(sprintf("%s: n=%d rho-GE vs r_g r=%.2f, %.0f%% sign-concordant\n",s,nrow(cc),r,100*mean(sign(cc[[paste0("rg_",s)]])==sign(cc[[paste0("ge_",s)]]))))}
cat("Panel A diseases:",uniqueN(rate$disease)," | r_g diseases:",nrow(rgt)," | rho-GE diseases:",nrow(rho),"\n")
print(rB[,.(method,sex,r,p=signif(p,2))])