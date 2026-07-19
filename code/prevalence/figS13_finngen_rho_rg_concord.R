## SUPP figure: FinnGen rho_GE (FusioS expression) vs r_g (LDSC genome-wide) concordance
## Shows the two disease--fertility correlation measures agree. Manuscript style, matches fig_finngen_rho.R.
suppressMessages({library(data.table);library(ggplot2);library(patchwork)})
FONT<-"Arial"; update_geom_defaults("text",list(family=FONT))
D<-"/Volumes/S840/04mypaper/lin/conflictresolution/rg_prevalence_analysis_2026-07-16/tables"
OUT<-c("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait",
       "/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/figures_final")
rho<-fread(file.path(D,"rhoGE_finngen_mhcx.tsv"))                 ## MHC-excluded rho_GE (full)
rg <-fread(file.path(D,"finngen_rg_vs_finregistry_prevalence.tsv"))[,.(phenocode,rg_father,rg_mother)]
d<-merge(rho,rg,by="phenocode")
m<-rbind(d[,.(phenocode,axis="Male",  rg=rg_father,rho=rho_father)],
         d[,.(phenocode,axis="Female",rg=rg_mother,rho=rho_mother)])[is.finite(rg)&is.finite(rho)]
m[,axisf:=factor(axis,levels=c("Male","Female"))]
BS<-8; COLp<-"#2C5985"; COLm<-"#C0392B"
th<-theme_classic(base_size=BS,base_family=FONT)+theme(text=element_text(family=FONT),
  axis.title=element_text(size=BS),axis.text=element_text(size=BS-1,color="grey20"),
  legend.text=element_text(size=BS-1.5),legend.title=element_blank(),legend.key.size=unit(3,"mm"),
  panel.grid=element_blank(),plot.tag=element_text(face="bold",size=11,family=FONT))
s<-m[,{c<-cor.test(rg,rho);.(r=c$estimate,p=c$p.value,n=.N)},by=axisf]; setorder(s,axisf)
lab<-s[,sprintf("%s: r=%.2f, P=%.0e (n=%d)",axisf,r,p,n)]
labcol<-s[,fifelse(axisf=="Male",COLp,COLm)]
print(s)
g<-ggplot(m,aes(rg,rho,color=axisf,fill=axisf))+
  geom_hline(yintercept=0,color="grey85",linetype=2,linewidth=.3)+geom_vline(xintercept=0,color="grey85",linetype=2,linewidth=.3)+
  geom_point(size=.6,alpha=.28)+
  geom_smooth(aes(linetype=axisf),method="lm",formula=y~x,se=TRUE,linewidth=.9,alpha=.15,color="black")+
  geom_smooth(aes(linetype=axisf),method="lm",formula=y~x,se=FALSE,linewidth=.8)+
  scale_color_manual(values=c(Male=COLp,Female=COLm),name=NULL)+
  scale_fill_manual(values=c(Male=COLp,Female=COLm),guide="none")+
  scale_linetype_manual(values=c(Male="solid",Female="22"),guide="none")+
  annotate("text",x=-Inf,y=Inf,hjust=-.05,vjust=c(1.4,2.9),label=lab,size=2.4,color=labcol,family=FONT)+
  labs(x=expression("disease--fertility "*italic(r)[g]*" (LDSC, genome-wide SNPs)"),
       y=expression("disease--fertility expression-effect ("*rho[GE]*")"))+
  th+theme(legend.position=c(.99,.02),legend.justification=c(1,0),
           legend.background=element_rect(fill=alpha("white",.7),color=NA))
for(o in OUT){ggsave(file.path(o,"fig_finngen_rho_rg_concord.pdf"),g,width=110,height=95,units="mm",device=cairo_pdf)
  ggsave(file.path(o,"fig_finngen_rho_rg_concord.png"),g,width=110,height=95,units="mm",dpi=300)}
cat(sprintf("\nfig_finngen_rho_rg_concord written (n=%d endpoints with rho+rg)\n",length(unique(m$phenocode))))
