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
P[,nrev:=round(revrate*ns)]

suppressMessages({library(ggplot2)});set.seed(1)
L<-list();ST<-list()
for(sc in c("Xf","Xm")){lab<-ifelse(sc=="Xf","Female","Male");D<-P[is.finite(get(sc))];D[,X:=get(sc)]
 D[,cls:=factor(fifelse(X<0,"opposite","same-sign"),levels=c("same-sign","opposite"))]
 disc<-D[X<0];conc<-D[X>0]
 a11<-sum(disc$nrev);a12<-sum(disc$ns)-a11;a21<-sum(conc$nrev);a22<-sum(conc$ns)-a21
 or<-(a11*a22)/(a12*a21);se<-sqrt(1/a11+1/a12+1/a21+1/a22);ci<-exp(log(or)+c(-1.96,1.96)*se)
 rawO<-100*a11/(a11+a12);rawS<-100*a21/(a21+a22)
 # disease-label permutation P on weighted opposite-effect difference
 obsdiff<-a11/(a11+a12)-a21/(a21+a22)
 rr<-if(sc=="Xf")rF else rM;gd<-rr[dz]
 perm<-replicate(2000,{gp<-setNames(sample(gd),dz);xp<-gp[D$a]*gp[D$b];i1<-xp<0;i2<-xp>0
   (sum(D$nrev[i1])/sum(D$ns[i1]))-(sum(D$nrev[i2])/sum(D$ns[i2]))})
 pp<-(1+sum(abs(perm)>=abs(obsdiff)))/2001
 D[,sex:=lab];L[[sc]]<-D
 pf<-function(p)ifelse(p<1e-3,sprintf("%.0e",p),sprintf("%.3f",p))
 ST[[sc]]<-data.table(sex=lab,box=sprintf("OR = %.2f [%.2f, %.2f]\nperm P = %s\nmean %.0f%% vs %.0f%%",or,ci[1],ci[2],pf(pp),rawO,rawS))}
A<-rbindlist(L);SB<-rbindlist(ST);mn<-A[,.(m=100*mean(revrate)),by=.(sex,cls)]
g<-ggplot(A,aes(cls,100*revrate,color=cls,fill=cls))+geom_violin(alpha=.13,color=NA,width=.9,scale="width")+
  geom_jitter(aes(size=ns),width=.11,alpha=.16,stroke=0)+geom_boxplot(width=.14,fill=NA,color="grey25",linewidth=.4,outlier.shape=NA)+
  geom_hline(yintercept=50,linetype=2,color="grey55",linewidth=.4)+geom_point(data=mn,aes(cls,m),inherit.aes=F,shape=23,size=3,fill="white",color="black",stroke=.8)+
  geom_text(data=SB,aes(x=1.5,y=116,label=box),inherit.aes=F,size=1.9,vjust=1,lineheight=.95)+
  facet_wrap(~sex,ncol=2)+scale_color_manual(values=c("same-sign"="#4575B4","opposite"="#D73027"),guide="none")+scale_fill_manual(values=c("same-sign"="#4575B4","opposite"="#D73027"),guide="none")+
  scale_size_continuous(range=c(.3,2.2),name="shared\ngenes")+scale_x_discrete(labels=function(v){m<-c(`same-sign`="atop(bold(\"fertility-concordant\"),rho[i]*rho[j]>0)",`opposite`="atop(bold(\"fertility-discordant\"),rho[i]*rho[j]<0)");parse(text=m[v])})+coord_cartesian(ylim=c(0,120))+
  labs(x=expression("fertility coupling of the two diseases   (sign of "*rho[list(GE,i)]*"\u00b7"*rho[list(GE,j)]*")"),y="shared genes with opposite\nexpression effect (%)")+
  theme_classic(base_size=7)+theme(strip.text=element_text(face="bold",size=7.5),strip.background=element_blank(),axis.text.x=element_text(size=6.2),axis.text.y=element_text(size=6),axis.title=element_text(size=7),axis.title.x=element_text(margin=margin(t=3)),legend.position="right",legend.key.size=unit(3,"mm"),legend.title=element_text(size=5.5),legend.text=element_text(size=5))
ggsave("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/Fig_intralocus_SUPP_box_sexmatched_withcancer.pdf",g,width=120,height=78,units="mm",device=cairo_pdf,bg="white")
ggsave("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/Fig_intralocus_SUPP_box_sexmatched_withcancer.png",g,width=120,height=78,units="mm",dpi=300,bg="white")
cat(sprintf("saved box. ndz=%d npairs=%d\n",length(dz),nrow(P)));print(SB)
