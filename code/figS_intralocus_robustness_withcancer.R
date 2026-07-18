suppressMessages({library(data.table);library(MASS);library(readxl)});options(width=200);set.seed(1)
strip<-function(x)sub("[.].*","",x);nk<-function(x)gsub("[^a-z0-9]","",tolower(x))
setwd("/Volumes/X10Pro/data/synergistic/mr_res")
st<-as.data.table(read_excel("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/Supplementary Table_v3.xlsx",sheet="Supplementary Table 1",skip=3));setnames(st,1:7,c("No","Trait","Type","Category","Source","Accession","PMID"));st<-st[!is.na(Trait)&grepl("^[0-9]",No)]
accmap<-setNames(st$Trait,gsub("[^A-Za-z0-9]","",st$Accession));nmap<-setNames(st$Trait,nk(st$Trait))
getacc<-function(x){m<-regmatches(x,regexpr("GCST[0-9]+",x));if(length(m))m else NA}
ov<-c(cad="Coronary Artery Disease",ibd="Inflammatory Bowel Disease",bip="Bipolar Disorder",adhd="ADHD",scz="Schizophrenia",ptsd="PTSD",Crohns="Crohn's Disease",Raynaud="Raynaud's Disease",Hyperplasiaofprostate="Benign Prostatic Hyperplasia",Highmyopia="High Myopia",osteoarthritis="Osteoarthritis",Heartfailure="Heart Failure",myocardialinfarction="Myocardial Infarction",majordepressivedisorder="Major Depressive Disorder",parkinson="Parkinson's Disease",alzheimer="Alzheimer's Disease",sle="Systemic Lupus Erythematosus",t1d="Type 1 Diabetes",atrialfibrillation="Atrial Fibrillation",endometriosis="Endometriosis",osteoporosis="Osteoporosis",hayfever="Hay Fever",bloodclot="Blood Clot",cardiometabolicmulti="Cardiometabolic Multimorbidity",hypertrophiccardiomyopathy="Hypertrophic Cardiomyopathy",metabolicsyndrome="Metabolic Syndrome",celiac="Celiac Disease",kidneystones="Kidney Stones")
idov<-c(Blood_Clot_Lung="Pulmonary Embolism",Chronic_obstructive_pulmonary_disease="Chronic Obstructive Pulmonary Disease")
mapname<-function(id){if(id%in%names(idov))return(unname(idov[id]));a<-getacc(id);if(!is.na(a)&&!is.na(accmap[a]))return(unname(accmap[a]));ck<-nk(gsub("GCST[0-9]+|EFO_?[0-9]+|\\.h\\.tsv\\.gz|\\.h\\.|_?k1$|[0-9]+$","",id));if(!is.na(nmap[ck]))return(unname(nmap[ck]));for(o in names(ov))if(grepl(nk(o),nk(id),fixed=T))return(ov[[o]]);id}
iscancer<-function(x)grepl("cancer|carcinoma|lymphoma|leukemia|myeloma|melanoma|Fibroid",x,ignore.case=T)
isContam<-function(x){z<-nk(x);z=="bloodclot"||grepl("celiac|sclerosis|ankylos|covid|narcoleps|metabolicsyndrome|cardiometabolic|manic|mentalhealth",z)}
co<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[,g:=strip(gene)];co<-co[chr%in%as.character(1:22)];chrm<-setNames(as.integer(co$chr),co$g);posm<-setNames(as.numeric(co$pos),co$g)
T5<-as.data.table(read_excel("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/Supplementary Table_v3.xlsx",sheet="Supplementary Table 5",skip=1));T5<-T5[!is.na(Disease)];T5[,k:=nk(Disease)]
## disease effects (sex-independent) from FEMALE ACAT file
d<-fread("ACAT_female_gene_disease_all.csv",select=c("gene","disease","b1_at_maxZ","padj_acat_1"))[,g:=strip(gene)][is.finite(b1_at_maxZ)]
d<-d[!(!is.na(chrm[g])&chrm[g]==6&posm[g]>25e6&posm[g]<34e6)]
DICT<-fread("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/disease_dictionary.csv")
d<-d[disease%in%DICT$acat_label];d[,ds:=sign(b1_at_maxZ)]
rf<-setNames(as.numeric(DICT$rhoGE_female),DICT$acat_label);rm<-setNames(as.numeric(DICT$rhoGE_male),DICT$acat_label)
d[,`:=`(rhoF=rf[disease],rhoM=rm[disease])]
dz<-unique(d[is.finite(rhoF)&is.finite(rhoM)]$disease);dz<-dz[!is.na(dz)]
rF<-setNames(unique(d[disease%in%dz,.(disease,rhoF)])$rhoF,unique(d[disease%in%dz,.(disease,rhoF)])$disease)
rM<-setNames(unique(d[disease%in%dz,.(disease,rhoM)])$rhoM,unique(d[disease%in%dz,.(disease,rhoM)])$disease)
## sign matrix (sig) + effect matrix (sim)
sig<-d[padj_acat_1<0.05&disease%in%dz,.(g,disease,ds)];Mx<-as.matrix(dcast(sig,g~disease,value.var="ds")[,-1]);rownames(Mx)<-dcast(sig,g~disease,value.var="ds")$g
EF<-as.matrix(dcast(d[disease%in%dz,.(g,disease,b1_at_maxZ)],g~disease,value.var="b1_at_maxZ")[,-1]);rownames(EF)<-dcast(d[disease%in%dz,.(g,disease,b1_at_maxZ)],g~disease,value.var="b1_at_maxZ")$g
sc<-colnames(Mx);out<-list()
for(i in 1:(length(dz)-1))for(j in (i+1):length(dz)){a<-dz[i];b<-dz[j];if(is.na(a)||is.na(b)||!(a%in%sc)||!(b%in%sc))next
  v<-Mx[,c(a,b)];v<-v[complete.cases(v),,drop=F];if(nrow(v)<15)next
  sim<-suppressWarnings(cor(EF[,a],EF[,b],use="pairwise.complete.obs"))
  out[[length(out)+1]]<-data.table(a,b,ns=nrow(v),revrate=mean(v[,1]!=v[,2]),sim=sim,Xf=rF[a]*rF[b],Xm=rM[a]*rM[b])}
P<-rbindlist(out)[is.finite(sim)];P[,Xmean:=(Xf+Xm)/2]
setnames(P,c('revrate','ns'),c('oep','n_shared'));P[,Xmn:=Xmean]
rmn<-setNames(sapply(dz,function(x)(rF[x]+rM[x])/2),dz)
library(ggplot2);library(patchwork)
cols<-c('Female (maternal)'='#C0392B','Male (paternal)'='#2C6FBB')

## ---- Panel A: pooled rank added-variable + disease-label permutation ----
Pa<-P[is.finite(Xmn)];Pa[,`:=`(rx=resid(lm(Xmn~sim,Pa)),ry=resid(lm(oep~sim,Pa)),Sm=sim)]
prA<-cor(Pa$rx,Pa$ry);perm<-replicate(10000,{gp<-setNames(sample(rmn[dz]),dz);cor(resid(lm(gp[Pa$a]*gp[Pa$b]~Pa$Sm)),Pa$ry)});pP<-(1+sum(abs(perm)>=abs(prA)))/10001
Pa[,dbin:=cut(rx,quantile(rx,seq(0,1,.1)),include.lowest=T,labels=1:10)];Bd<-Pa[,.(my=mean(ry),se=sd(ry)/sqrt(.N),n=.N),by=dbin][!is.na(dbin)][order(dbin)]
pA<-ggplot(Bd,aes(as.integer(dbin),my))+geom_hline(yintercept=0,color="grey85",linewidth=.3)+
  geom_smooth(method="lm",color="#C0392B",fill="#C0392B",alpha=.12,linewidth=.6)+geom_errorbar(aes(ymin=my-se,ymax=my+se),width=0,linewidth=.5,color="grey45")+geom_point(size=2.1,color="#1A5276")+annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=1.4,label=sprintf("n = %d pairs / bin",round(nrow(Pa)/10)),size=2.2,color="grey45")+scale_x_continuous(breaks=c(1,5.5,10),labels=c("fertility-\ndiscordant","decile","fertility-\nconcordant"),expand=expansion(mult=.09))+
  labs(x="how the two diseases couple to fertility",y="opposite-effect proportion\n(residual, mean \u00b1 SE)")+
    annotate("text",x=-Inf,y=-Inf,hjust=-0.06,vjust=-0.55,label=sprintf("partial r = %.2f\npermutation P = %.4f\n(disease sim. controlled)",prA,pP),size=2.3,color="grey30",lineheight=.95)+
  theme_classic(8)+theme(axis.title=element_text(size=8),axis.text.x=element_text(size=6.3),plot.margin=margin(4,8,4,4))+coord_cartesian(clip="off")
## ---- Panels B-D: sex-specific standardized beta robustness ----
sl<-function(d){m<-lm(Y~X+sim,d);c(coef(m)["X"],confint(m)["X",],summary(m)$coef["X","Pr(>|t|)"])}
res_lk<-list();res_lodo<-list();ROB<-list()
for(sc in c("Xf","Xm")){sx<-ifelse(sc=="Xf","Female (maternal)","Male (paternal)");D<-P[is.finite(get(sc))];D[,`:=`(X=scale(get(sc))[,1],Y=scale(oep)[,1])];D[,cook:=cooks.distance(lm(Y~X+sim,D))]
 for(k in c(0,1,5,10,20,50)){dd<-if(k==0)D else D[order(-cook)][-(1:k)];sk<-sl(dd);res_lk[[length(res_lk)+1]]<-data.table(sex=sx,k=k,slope=sk[1],lo=sk[2],hi=sk[3],P=sk[4])}
 hub<-summary(rlm(Y~X+sim,D,maxit=200))$coef["X",];mmm<-summary(rlm(Y~X+sim,D,method="MM",maxit=200))$coef["X",]
 q<-quantile(D$cook,.99);mt<-lm(Y~X+sim,D[cook<=q]);s0<-sl(D);ct<-confint(mt)["X",]
 ROB[[length(ROB)+1]]<-data.table(sex=sx,est=c("OLS","Huber","MM","1%-trim"),slope=c(s0[1],hub[1],mmm[1],coef(mt)["X"]),lo=c(s0[2],hub[1]-1.96*hub[2],mmm[1]-1.96*mmm[2],ct[1]),hi=c(s0[3],hub[1]+1.96*hub[2],mmm[1]+1.96*mmm[2],ct[2]))
 for(d in dz){dd<-D[a!=d&b!=d];if(nrow(dd)<50)next;res_lodo[[length(res_lodo)+1]]<-data.table(sex=sx,slope=coef(lm(Y~X+sim,dd))["X"])}}
LK<-rbindlist(res_lk);ROB<-rbindlist(ROB);LODO<-rbindlist(res_lodo);full<-LK[k==0]
data.table::fwrite(LK,"/Volumes/S840/04mypaper/lin/conflictresolution/fig4_intralocus_session_20260714/TableS_leaveTopK_cooks.csv")
cat("=ROB=\n");print(ROB[,.(sex,est,slope=round(slope,3),lo=round(lo,3),hi=round(hi,3))]);cat("=LODO=\n");print(LODO[,.(n=.N,neg=sum(slope<0),med=round(median(slope),3),mn=round(min(slope),3),mx=round(max(slope),3)),by=sex]);cat("=FULL=\n");print(full[,.(sex,beta=round(slope,3))])
pB<-ggplot(LK,aes(factor(k),slope,color=sex,group=sex))+geom_hline(yintercept=0,linetype=3,color="grey60")+geom_hline(data=full,aes(yintercept=slope,color=sex),linetype=2,linewidth=.5,show.legend=F)+
  geom_line(linewidth=.7)+geom_point(size=2)+scale_color_manual(values=cols,name=NULL)+labs(x="pairs removed (Cook's D)",y="partial regression coefficient (β)")+
  theme_classic(8)+theme(plot.title=element_text(size=8.8),axis.title=element_text(size=8))
ROB[,est:=factor(est,levels=c("1%-trim","MM","Huber","OLS"))]
pC<-ggplot(ROB,aes(slope,est,color=sex))+geom_vline(xintercept=0,linetype=3,color="grey60")+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.22,linewidth=.6,position=position_dodge(.55))+
  geom_point(size=2,position=position_dodge(.55))+scale_color_manual(values=cols,name=NULL)+guides(color="none")+labs(x="partial regression coefficient (β)  [95% CI]",y=NULL)+
  theme_classic(8)+theme(plot.title=element_text(size=8.8),axis.title=element_text(size=8))
pD<-ggplot(LODO,aes(slope,fill=sex))+geom_vline(xintercept=0,linetype=3,color="grey60")+geom_histogram(bins=22,alpha=.6,color="white",linewidth=.12,position="identity")+
  geom_vline(data=full,aes(xintercept=slope,color=sex),linetype=2,linewidth=.6,show.legend=F)+scale_fill_manual(values=cols,name=NULL)+scale_color_manual(values=cols,name=NULL)+guides(fill="none",color="none")+
  labs(x="partial regression coefficient\n(β, each disease removed once)",y="diseases (count)")+theme_classic(8)+theme(plot.title=element_text(size=8.8),axis.title=element_text(size=8))
fig<-(pA+pC+pD+plot_layout(nrow=1,guides="collect"))+plot_annotation(tag_levels="A") & theme(legend.position="bottom",legend.justification="center")
BD<-"/Volumes/S840/04mypaper/lin/conflictresolution/fig4_intralocus_session_20260714/"
ggsave(paste0(BD,"FigS_intralocus_robustness_3panel_withcancer.png"),fig,width=245,height=88,units="mm",dpi=300,bg="white")
ggsave(paste0(BD,"FigS_intralocus_robustness_3panel_withcancer.pdf"),fig,width=245,height=88,units="mm",device=cairo_pdf,bg="white")
cat(sprintf("saved 4-panel. partialSpearman=%.3f permP=%.4f\n",prA,pP))
