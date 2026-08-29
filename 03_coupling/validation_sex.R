suppressMessages({library(data.table);library(ggplot2);library(patchwork)})
MR<-"${MR}"; strip<-function(x)sub("[.].*","",x); nk<-function(x)gsub("[^a-z0-9]","",tolower(x))
co<-fread(file.path(MR,"gene_coords_v19.tsv"),col.names=c("g","chr","pos"))[,g:=strip(g)]
rec<-rbindlist(lapply(1:22,function(c){r<-fread(file.path(MR,"recomb",sprintf("plink.chr%d.GRCh37.map",c)),header=F);data.table(chr=c,bp=r$V4,cm=r$V3)}))
co[,cm:=NA_real_];for(c in 1:22){i<-which(co$chr==c);m<-rec[chr==c];if(length(i))co$cm[i]<-approx(m$bp,m$cm,co$pos[i],rule=2,ties="ordered")$y}
isM<-function(chr,pos)chr==6&pos>=25e6&pos<=34e6
b1<-fread(file.path(MR,"gene_level_b1_disease_acat_summary.dropbox.csv"),select=c("gene","disease","top_beta","top_se"))[,g:=strip(gene)][is.finite(top_beta)&top_se>0][,z1:=top_beta/top_se];b1<-merge(b1,co,by="g")[!isM(chr,pos)&abs(z1)<9]
b2<-fread(file.path(MR,"gene_level_b2_fitness_acat_summary.dropbox.csv"),select=c("gene","fitness","top_beta","top_se"))[,g:=strip(gene)][is.finite(top_beta)&top_se>0][,z2:=top_beta/top_se];b2<-merge(b2,co,by="g")[!isM(chr,pos)&abs(z2)<9]
rhoAll<-function(sx){f<-b2[fitness==sx];rbindlist(lapply(unique(b1$disease),function(dz){m<-merge(b1[disease==dz,.(g,z1,chr,cm)],f[,.(g,z2)],by="g");if(nrow(m)<30)return(NULL);bk<-m[order(-abs(z1)),.SD[1],by=paste(chr,floor(cm/1))];data.table(dk=nk(dz),all=cor(m$z1,m$z2),blk=cor(bk$z1,bk$z2))}))}
R<-rbind(rhoAll("male")[,sex:="Male"],rhoAll("female")[,sex:="Female"])
# null per sex
set.seed(1);nl<-rbindlist(lapply(c("male","female"),function(sx){f<-b2[fitness==sx];v<-c();for(i in 1:15){f2<-copy(f)[,z2:=sample(z2)];for(dz in sample(unique(b1$disease),20)){m<-merge(b1[disease==dz,.(g,z1)],f2[,.(g,z2)],by="g");if(nrow(m)>=30)v<-c(v,cor(m$z1,m$z2))}};data.table(rho=v,sex=fifelse(sx=="male","Male","Female"))}))
# LDSC per sex
Lr<-readLines(file.path(MR,"session_2026-07-09_iv_filter","dz_rg_97_ALLESTIMABLE.txt"))
vv<-function(p,t){mm<-grep(paste0("^",t,"="),p,value=T);if(!length(mm))return(NA);suppressWarnings(as.numeric(sub(paste0("^",t,"=[[:space:]]*(-?[0-9.eE-]+).*"),"\\1",mm)))}
rg<-rbindlist(lapply(Lr,function(x){pp<-strsplit(x,"\t")[[1]];data.table(dk=nk(sub("_$","",pp[1])),Male=vv(pp,"father"),Female=vv(pp,"mother"))}))
rgl<-melt(rg,id.vars="dk",variable.name="sex",value.name="ldsc")[!is.na(ldsc)][,sex:=as.character(sex)]
R<-merge(R,rgl,by=c("dk","sex"),all.x=T)
## Table S5 covers 90 diseases; merge_tableS5_90.R builds it from the two partial
## subsets (rhoge_auth_compare.csv, 61; the earlier 77-disease table).
comp<-fread("${PROJ}/tableS5_90diseases.tsv")
compS<-rbind(comp[,.(dk=Disease,mr=rhoGE_FusioMR_male,fs=rhoGE_FusioS_male,sex="Male")],
             comp[,.(dk=Disease,mr=rhoGE_FusioMR_female,fs=rhoGE_FusioS_female,sex="Female")])
cols<-c("Male"="#2C5985","Female"="#C0392B")
lb<-function(d,x,y)d[,{r<-cor(get(x),get(y),use="complete");sprintf("%s: r=%.2f",sex,r)},by=sex]
th<-theme_classic(base_size=9)+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7.5,color="grey35"),legend.position="top",legend.title=element_blank())
scatter<-function(d,x,y,xl,yl,ti,diag=F){l<-setNames(lb(d[is.finite(get(x))&is.finite(get(y))],x,y)$V1,lb(d[is.finite(get(x))&is.finite(get(y))],x,y)$sex)
  p<-ggplot(d,aes(.data[[x]],.data[[y]],color=sex))+geom_hline(yintercept=0,linetype=2,color="grey80")+geom_vline(xintercept=0,linetype=2,color="grey80")
  if(diag)p<-p+geom_abline(linetype=3,color="grey60")
  p+geom_point(size=1.5,alpha=.7)+geom_smooth(method="lm",se=F,linewidth=.5)+scale_color_manual(values=cols,labels=l)+labs(x=xl,y=yl,title=ti)+th}
## A (moved from main Fig): genome-wide SNP r_g ALSO predicts % antagonistic loci, by sex
CF<-fread(file.path(MR,"session_2026-07-09_iv_filter","colleague_shared_iv_filter.csv"))[keep_after_filter==TRUE]
rate<-CF[,.(n=.N,pct=100*mean(combo%in%c("++","--"))),by=.(dk=nk(disease),sex=fifelse(fitness=="male","Male","Female"))][n>=1]
Aant<-merge(rate,rgl,by=c("dk","sex"))[is.finite(ldsc)]
rA<-Aant[,.(lab=sprintf("%s: r=%.2f",sex,cov.wt(cbind(ldsc,pct),wt=n,cor=T)$cor[1,2])),by=sex]
pA<-ggplot(Aant,aes(ldsc,pct,color=sex))+geom_vline(xintercept=0,linetype=2,color="grey80")+geom_point(aes(size=n),alpha=.6)+geom_smooth(aes(weight=n),method="lm",formula=y~x,se=T,linewidth=.6)+
  scale_size_continuous(range=c(.5,2.6),name="indep. loci",breaks=c(1,5,20,50))+scale_color_manual(values=cols)+
  geom_text(data=rA[sex=="Male"],aes(x=-Inf,y=Inf,label=lab),color="#2C5985",hjust=-.08,vjust=1.5,size=2.6,inherit.aes=F)+geom_text(data=rA[sex=="Female"],aes(x=-Inf,y=Inf,label=lab),color="#C0392B",hjust=-.08,vjust=3,size=2.6,inherit.aes=F)+
  labs(x=expression(italic(r)[g]*" (LDSC, genome-wide SNPs)"),y="% antagonistic (indep. loci)")+th
pB<-scatter(compS,"mr","fs",expression("FusioMR "*rho[GE]*" (shared IV)"),expression("FusioS "*rho[GE]*" (independent IV)"),NULL,diag=T)
pC<-scatter(R,"all","blk",expression(rho[GE]*" (all genes)"),expression(rho[GE]*" (cM-block)"),NULL,diag=T)
OUT<-"${PROJ}/figures_final"
## Figure 1 (standalone): SNP r_g also predicts antagonism
ggsave(file.path(OUT,"FigS_rg_antagonism.pdf"),pA,width=5.4,height=4.6,device=cairo_pdf)
ggsave(file.path(OUT,"FigS_rg_antagonism.png"),pA,width=5.4,height=4.6,dpi=160)
## Figure 2 (A|B): rho_GE robust to between-trait & between-gene shared IVs
figR<-(pB|pC)+plot_annotation(tag_levels="a")&theme(plot.tag=element_text(face="bold",size=13))
ggsave(file.path(OUT,"FigS_rhoGE_robustness.pdf"),figR,width=10,height=4.6,device=cairo_pdf)
ggsave(file.path(OUT,"FigS_rhoGE_robustness.png"),figR,width=10,height=4.6,dpi=160)
cat("saved 2 figs. r values:\n");print(lb(Aant,"ldsc","pct"));print(lb(compS,"mr","fs"));print(lb(R,"all","blk"))
