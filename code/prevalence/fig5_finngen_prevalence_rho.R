## FinnGen rho_GE (FusioS expression-effect) vs FinRegistry prevalence — manuscript style
## MAIN-text figure (rho version). Supplement keeps the rg version (fig_finngen_prevalence.pdf).
suppressMessages({library(data.table);library(ggplot2);library(patchwork)})
FONT<-"Arial"; update_geom_defaults("text",list(family=FONT))
D<-"/Volumes/S840/04mypaper/lin/conflictresolution/rg_prevalence_analysis_2026-07-16/tables"
OUT<-c("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait",
       "/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/figures_final")
## rho_GE per FinnGen disease (father/mother) — freshest pull from Randi
rho<-fread(file.path(D,"rhoGE_finngen_mhcx.tsv"))  ## MHC-excluded (chr6:25-34Mb)
## FinRegistry prevalence (log10 % ) keyed by phenocode
prev<-fread(file.path(D,"finngen_rg_vs_finregistry_prevalence.tsv"))[,.(phenocode,lp=log10_finregistry_prev_pct)]
d<-merge(rho,prev,by="phenocode")
m<-melt(d[,.(phenocode,lp,Male=rho_father,Female=rho_mother)],id.vars=c("phenocode","lp"),
        variable.name="axis",value.name="rho")[is.finite(rho)&is.finite(lp)]
m[,axisf:=factor(axis,levels=c("Male","Female"))]
BS<-8; COLp<-"#2C5985"; COLm<-"#C0392B"
th<-theme_classic(base_size=BS,base_family=FONT)+theme(text=element_text(family=FONT),
  axis.title=element_text(size=BS),axis.text=element_text(size=BS-1,color="grey20"),
  legend.text=element_text(size=BS-1.5),legend.title=element_blank(),legend.key.size=unit(3,"mm"),
  panel.grid=element_blank(),plot.tag=element_text(face="bold",size=11,family=FONT))
s<-m[,{c<-cor.test(rho,lp);.(r=c$estimate,p=c$p.value,n=.N)},by=axisf]
setorder(s,axisf)
lab<-s[,sprintf("%s: r=%.2f, P=%.0e (n=%d)",axisf,r,p,n)]
labcol<-s[,fifelse(axisf=="Male",COLp,COLm)]
print(s)
## Panel A: rho_GE vs prevalence, paternal & maternal
pA<-ggplot(m,aes(rho,lp,color=axisf,fill=axisf))+geom_vline(xintercept=0,color="grey85",linetype=2,linewidth=.3)+
  geom_point(size=.6,alpha=.30)+
  geom_smooth(aes(linetype=axisf),method="lm",formula=y~x,se=TRUE,linewidth=.9,alpha=.15,color="black")+
  geom_smooth(aes(linetype=axisf),method="lm",formula=y~x,se=FALSE,linewidth=.8)+
  scale_color_manual(values=c(Male=COLp,Female=COLm),name=NULL)+
  scale_fill_manual(values=c(Male=COLp,Female=COLm),guide="none")+
  scale_linetype_manual(values=c(Male="solid",Female="22"),guide="none")+
  annotate("text",x=-Inf,y=Inf,hjust=-.05,vjust=c(1.4,2.9),label=lab,size=2.2,color=labcol,family=FONT)+
  labs(x=expression("disease--fertility expression-effect ("*rho[GE]*")"),
       y=expression(log[10]*" FinRegistry prevalence (%)"))+
  th+theme(legend.position=c(.99,.02),legend.justification=c(1,0),legend.background=element_rect(fill=alpha("white",.7),color=NA))
## Panel B: rho_GE>0 vs <0 prevalence, by fertility axis
m[,grp:=factor(fifelse(rho>0,"rho > 0","rho < 0"),levels=c("rho < 0","rho > 0"))]
st<-m[,{w<-wilcox.test(lp~grp);.(p=w$p.value,y=max(lp)+.3)},by=axisf]; print(st)
pB<-ggplot(m,aes(grp,lp))+geom_violin(aes(fill=grp),color=NA,alpha=.24,scale="width",width=.85)+
  geom_boxplot(width=.16,outlier.size=.25,outlier.alpha=.15,linewidth=.35,fill="white")+
  facet_wrap(~axisf)+scale_fill_manual(values=c(`rho < 0`=COLp,`rho > 0`=COLm),guide="none")+
  geom_text(data=st,aes(x=1.5,y=y,label=sprintf("italic(P)=='%.0e'",p)),parse=TRUE,size=2.3,inherit.aes=FALSE,family=FONT)+
  labs(x=expression("disease--fertility "*rho[GE]),y=expression(log[10]*" FinRegistry prevalence (%)"))+
  theme_bw(BS,base_family=FONT)+theme(text=element_text(family=FONT),panel.grid=element_blank(),
    strip.background=element_rect(fill="grey92",color=NA),strip.text=element_text(face="bold",size=BS-1),
    axis.text=element_text(size=BS-1,color="grey20"),plot.tag=element_text(face="bold",size=11,family=FONT))
fig<-pA/pB+plot_layout(heights=c(1.05,1))+plot_annotation(tag_levels=list(c("A","B")))&theme(plot.tag=element_text(face="bold",size=11,family=FONT))
for(o in OUT){ggsave(file.path(o,"fig_finngen_prevalence_rho.pdf"),fig,width=120,height=165,units="mm",device=cairo_pdf)
  ggsave(file.path(o,"fig_finngen_prevalence_rho.png"),fig,width=120,height=165,units="mm",dpi=300)}
cat(sprintf("\nfig_finngen_prevalence_rho written (n endpoints with rho+prev = %d)\n",length(unique(m$phenocode))))
