# NEW Fig 4 = FIVE panels (all native, canonical repo code; numbers unchanged):
#  a: lollipop (MALE) | b: rhoGE-vs-rg concordance (both sexes) | c: %antag-vs-rhoGE (both sexes, NO band) | d: triangle (MALE) | e: violin (MALE). female d/e -> supplement
pmP<-function(p) ifelse(p<1e-3, sprintf("%.1f%%*%%10^%d", p/10^floor(log10(p)), floor(log10(p))), sprintf("'%.2g'", p))
suppressMessages({library(data.table);library(ggplot2);library(patchwork);library(cowplot);library(readxl)});options(width=200);set.seed(1)
setwd("${MR}")
CFPATH<-"${MR}/session_2026-07-09_iv_filter/colleague_shared_iv_filter.csv"
MY<-"${ROOT}"
STv3<-"${PROJ}/Supplementary Table_v3.xlsx"
DICTP<-"${PROJ}/disease_dictionary.csv"
BS<-8; th<-theme_classic(base_size=BS)+theme(axis.title=element_text(size=BS),axis.text=element_text(size=BS-1,color="grey20"),legend.text=element_text(size=BS-1),legend.key.size=unit(3.2,"mm"),panel.grid=element_blank())
nk<-function(x)gsub("[^a-z0-9]","",tolower(x)); strip<-function(x)sub("[.].*","",x)
co<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[chr%in%as.character(1:22)][,.(gene,chr=as.integer(chr),pos=as.numeric(pos))]
isMHC<-function(chr,pos) chr==6 & pos>25e6 & pos<34e6
## ---- data (fig3_plot_panels.R) ----
b1<-fread("gene_level_b1_disease_acat_summary.dropbox.csv",select=c("gene","disease","top_beta","top_se","p_acat"))
b2<-fread("gene_level_b2_fitness_acat_summary.dropbox.csv",select=c("gene","fitness","top_beta","top_se","p_acat"))[!is.na(fitness)]
za<-merge(b1[is.finite(top_beta)&top_se>0,.(gene,disease,z1=top_beta/top_se)],co,by="gene")[!isMHC(chr,pos)&abs(z1)<9]
zb<-merge(b2[is.finite(top_beta)&top_se>0,.(gene,sex=fifelse(fitness=="male","M","F"),z2=top_beta/top_se)],co,by="gene")[!isMHC(chr,pos)&abs(z2)<9]
rge<-merge(za[,.(gene,disease,z1)],zb[,.(gene,sex,z2)],by="gene",allow.cartesian=T)[,{ct<-suppressWarnings(cor.test(z1,z2));.(rge=unname(ct$estimate),ng=.N,p_ge=ct$p.value)},by=.(disease,sex)][ng>=30]
rho<-dcast(rge,disease~sex,value.var="rge"); rhop<-dcast(rge,disease~sex,value.var="p_ge");setnames(rhop,c("M","F"),c("pgeM","pgeF"))
Lr<-readLines("dz_rg_97_ALLESTIMABLE.txt");vv<-function(p,t){m<-grep(paste0("^",t,"="),p,value=T);if(!length(m))return(c(NA,NA));c(suppressWarnings(as.numeric(sub(paste0("^",t,"=[[:space:]]*(-?[0-9.eE-]+).*"),"\\1",m))),suppressWarnings(as.numeric(sub(".*\\(([0-9.eE-]+)\\).*","\\1",m))))}
rgt<-rbindlist(lapply(Lr,function(x){p<-strsplit(x,"\t")[[1]];fa<-vv(p,"father");mo<-vv(p,"mother");data.table(disease=sub("_$","",p[1]),M=fa[1],se_fa=fa[2],F=mo[1],se_mo=mo[2])}))[!is.na(disease)&disease!=""&is.finite(M)&is.finite(F)][,`:=`(zM=M/se_fa,zF=F/se_mo)]
CF<-fread(CFPATH)[keep_after_filter==TRUE]
rate<-CF[,.(n=.N,pct=100*mean(combo%in%c("++","--"))),by=.(dk=nk(disease),sex=fifelse(fitness=="male","M","F"))][n>=1]
st1<-as.data.table(read_excel("Supplementary_Table_1.xlsx",sheet=1,skip=2));setnames(st1,2:4,c("Trait","Type","Category"));st1<-st1[grepl("Disease",Type,ignore.case=T),.(k=nk(Trait),Category)]
brc<-fread("Supplementary_Table_1_corrected_v2.csv",select=c("b1_disease_id","Suppl_Trait"));brc[,k:=nk(Suppl_Trait)];brc<-merge(brc,st1,by="k",all.x=T)
map6<-function(c)fifelse(c=="Neuropsychiatric","Neuropsychiatric",fifelse(grepl("Immune",c),"Immune",fifelse(c%in%c("Cardiovascular","Metabolic"),"Cardiometabolic",fifelse(c=="Cancers","Cancer",fifelse(c=="Reproductive","Reprod","Other")))))
brc[,cat6:=map6(Category)]
CO<-rbind(rge[,.(dk=nk(disease),sex,corr=rge,method="expression effect (rho-GE, MR)")],melt(rgt[,.(dk=nk(disease),M,F)],id.vars="dk",variable.name="sex",value.name="corr")[is.finite(corr)][,`:=`(sex=as.character(sex),method="genome-wide SNPs (r_g, LDSC)")])
db<-merge(rate,CO,by=c("dk","sex"),allow.cartesian=T);db[,method:=factor(method,levels=c("expression effect (rho-GE, MR)","genome-wide SNPs (r_g, LDSC)"))];db[,sex:=factor(sex,levels=c("M","F"),labels=c("Male","Female"))]
rB<-db[,{w<-n;rw<-suppressWarnings(cov.wt(cbind(corr,pct),wt=w,cor=TRUE)$cor[1,2]);pw<-tryCatch(coef(summary(lm(pct~corr,weights=w)))[2,4],error=function(e)NA);.(r=round(rw,2),p=pw)},by=.(method,sex)];rB[,lab:=sprintf("'%s: r=%.2f, '*italic(P)==%s",ifelse(sex=="Male","M","F"),r,pmP(p))]
dbL<-db[method=="expression effect (rho-GE, MR)"];rBL<-rB[method=="expression effect (rho-GE, MR)"]
COLS<-c(Male="#2C5985",Female="#C0392B"); LTS<-c(Male="solid",Female="22")
cat("PANEL B size variable n: median", median(dbL$n), " range", range(dbL$n), "\n")
## the rendered radius under scale_size_continuous(range=c(.35,2)) for the median point
cat("  quartiles of n:", quantile(dbL$n, c(.25,.5,.75)), "\n")
sz <- .35 + (2-.35)*sqrt((dbL$n-min(dbL$n))/(max(dbL$n)-min(dbL$n)))
cat("  rendered sizes: median", round(median(sz),2), " quartiles", round(quantile(sz,c(.25,.75)),2), "\n")
cat("  -> equivalent fixed size:", round(.35 + (2-.35)*sqrt((median(dbL$n)-min(dbL$n))/(max(dbL$n)-min(dbL$n))), 2), "\n")

## ---- panel c = %antag vs rhoGE, NO band ----
pc<-ggplot(dbL,aes(corr,pct,color=sex,fill=sex))+geom_point(aes(size=n),alpha=.5)+geom_smooth(aes(weight=n,linetype=sex),method="lm",formula=y~x,se=FALSE,linewidth=.8)+
  geom_text(data=rBL[sex=="Male"],aes(x=-Inf,y=Inf,label=lab),color=COLS["Male"],hjust=-.06,vjust=1.5,size=2.3,parse=TRUE,inherit.aes=F)+geom_text(data=rBL[sex=="Female"],aes(x=-Inf,y=Inf,label=lab),color=COLS["Female"],hjust=-.06,vjust=3,size=2.3,parse=TRUE,inherit.aes=F)+
  scale_size_continuous(range=c(.35,2),name="indep. loci",breaks=c(1,5,20,50))+scale_color_manual(values=COLS,name=NULL)+scale_fill_manual(values=COLS,name=NULL)+scale_linetype_manual(values=LTS,name=NULL)+guides(fill="none")+
  labs(x=expression("disease-fertility "*rho[GE]*" (expression effect)"),y="Antagonistic gene fraction (%, indep. loci)")+th+theme(legend.position=c(.995,.02),legend.justification=c(1,0),legend.background=element_rect(fill=alpha("white",.7),color=NA),legend.box="horizontal")
## ---- panel b = rhoGE vs rg concordance ----
conv<-merge(rgt[,.(dk=nk(disease),rgM=M,rgF=F)],rho[,.(dk=nk(disease),geM=M,geF=F)],by="dk")
convL<-rbind(conv[,.(sex="Male",rg=rgM,ge=geM)],conv[,.(sex="Female",rg=rgF,ge=geF)])[is.finite(rg)&is.finite(ge)][,sex:=factor(sex,levels=c("Male","Female"))]
rC<-convL[,{ct<-suppressWarnings(cor.test(rg,ge));.(lab=sprintf("'%s: r=%.2f, '*italic(P)==%s",sex,ct$estimate,pmP(ct$p.value)))},by=sex]
pb<-ggplot(convL,aes(rg,ge,color=sex))+geom_hline(yintercept=0,linetype=2,color="grey85")+geom_vline(xintercept=0,linetype=2,color="grey85")+
  geom_point(size=1,alpha=.6)+geom_smooth(aes(linetype=sex),method="lm",formula=y~x,se=F,linewidth=.7)+
  geom_text(data=rC[sex=="Male"],aes(x=-Inf,y=Inf,label=lab),color=COLS["Male"],hjust=-.08,vjust=1.5,size=2.3,parse=TRUE,inherit.aes=F)+geom_text(data=rC[sex=="Female"],aes(x=-Inf,y=Inf,label=lab),color=COLS["Female"],hjust=-.08,vjust=3,size=2.3,parse=TRUE,inherit.aes=F)+
  scale_color_manual(values=COLS,name=NULL)+scale_linetype_manual(values=LTS,name=NULL)+labs(x=expression(italic(r)[g]*" (LDSC, genome-wide SNPs)"),y=expression(rho[GE]*" (MR, expression effect)"))+th+theme(legend.position="none")
## ---- panel a = lollipop, MALE only ----
aC<-merge(rgt[,.(disease,rg_M=M,rg_F=F,zM,zF)],rho[,.(disease,ge_M=M,ge_F=F)],by="disease",all.x=T);aC<-merge(aC,rhop[,.(disease,pgeM,pgeF)],by="disease",all.x=T)
aC[,k:=nk(disease)];aC<-merge(aC,brc[,.(k2=nk(b1_disease_id),cat6,Suppl_Trait)],by.x="k",by.y="k2",all.x=T)
aC[is.na(cat6),cat6:="Other"];aC[is.na(Suppl_Trait),Suppl_Trait:=disease];aC[,Suppl_Trait:=gsub("_"," ",Suppl_Trait)];aC[,ord:=rowMeans(cbind(rg_M,rg_F),na.rm=T)]
rgcat<-aC[,.(m=mean(ord,na.rm=T)),by=cat6][order(-m)];aC[,cat6:=factor(cat6,levels=c(setdiff(rgcat$cat6,"Other"),"Other"))]
aC<-aC[order(cat6,ord)];aC[,xpos:=1:.N];bd<-aC[,.(mid=mean(xpos),x1=max(xpos)),by=cat6][order(cat6)];seps<-head(bd$x1,-1)+.5
dCm<-rbindlist(list(aC[,.(Suppl_Trait,cat6,xpos,measure="r_g",corr=rg_M,sig=abs(zM)>1.96)],aC[,.(Suppl_Trait,cat6,xpos,measure="rho_GE",corr=ge_M,sig=pgeM<0.05)]))
dCm[,corr_s:=corr/max(abs(corr),na.rm=T),by=measure];dCm[,xplot:=xpos+ifelse(measure=="r_g",-.22,.22)]
dCm[,grp:=fifelse(!sig|is.na(sig),"n.s.",fifelse(measure=="r_g","genome-wide SNPs (r_g, LDSC)","expression effect (rho-GE)"))]
dCm[,grp:=factor(grp,levels=c("genome-wide SNPs (r_g, LDSC)","expression effect (rho-GE)","n.s."))]
cols<-c("genome-wide SNPs (r_g, LDSC)"="#2C5985","expression effect (rho-GE)"="#C0392B","n.s."="grey65");shp<-c("genome-wide SNPs (r_g, LDSC)"=16,"expression effect (rho-GE)"=16,"n.s."=1)
pa<-ggplot(dCm,aes(xplot,corr_s,color=grp,shape=grp))+geom_hline(yintercept=0,color="grey55",linewidth=.3)+geom_vline(xintercept=seps,linetype="dashed",color="grey70",linewidth=.3)+geom_linerange(aes(ymin=0,ymax=corr_s),linewidth=.3)+geom_point(size=.9,stroke=.4)+
  geom_text(data=bd,aes(x=mid,y=.92,label=cat6),inherit.aes=FALSE,fontface="bold",size=2)+scale_color_manual(values=cols,name=NULL)+scale_shape_manual(values=shp,name=NULL)+
  scale_x_continuous(breaks=aC$xpos,labels=aC$Suppl_Trait,expand=expansion(add=.6))+coord_cartesian(ylim=c(-1.03,1.03),clip="off")+labs(x=NULL,y="disease-fertility correlation\n(scaled, male axis)")+
  theme_bw(7)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=4.2),panel.grid=element_blank(),legend.position="top",legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=5.6),legend.margin=margin(0,0,4,0),plot.margin=margin(1,4,10,1))

## ================= panels d (triangle) + e (violin): fig4_AB_global3bin.R =================
st<-as.data.table(read_excel(STv3,sheet="Supplementary Table 1",skip=3));setnames(st,1:7,c("No","Trait","Type","Category","Source","Accession","PMID"));st<-st[!is.na(Trait)&grepl("^[0-9]",No)]
accmap<-setNames(st$Trait,gsub("[^A-Za-z0-9]","",st$Accession));nmap<-setNames(st$Trait,nk(st$Trait))
getacc<-function(x){m<-regmatches(x,regexpr("GCST[0-9]+",x));if(length(m))m else NA}
ov<-c(cad="Coronary Artery Disease",ibd="Inflammatory Bowel Disease",bip="Bipolar Disorder",adhd="ADHD",scz="Schizophrenia",ptsd="PTSD",Crohns="Crohn's Disease",Raynaud="Raynaud's Disease",Hyperplasiaofprostate="Benign Prostatic Hyperplasia",Highmyopia="High Myopia",osteoarthritis="Osteoarthritis",Heartfailure="Heart Failure",myocardialinfarction="Myocardial Infarction",majordepressivedisorder="Major Depressive Disorder",parkinson="Parkinson's Disease",alzheimer="Alzheimer's Disease",sle="Systemic Lupus Erythematosus",t1d="Type 1 Diabetes",atrialfibrillation="Atrial Fibrillation",endometriosis="Endometriosis",osteoporosis="Osteoporosis",hayfever="Hay Fever",bloodclot="Blood Clot",cardiometabolicmulti="Cardiometabolic Multimorbidity",hypertrophiccardiomyopathy="Hypertrophic Cardiomyopathy",metabolicsyndrome="Metabolic Syndrome",celiac="Celiac Disease",kidneystones="Kidney Stones")
idov<-c(Blood_Clot_Lung="Pulmonary Embolism",Chronic_obstructive_pulmonary_disease="Chronic Obstructive Pulmonary Disease")
mapname<-function(id){if(id%in%names(idov))return(unname(idov[id]));a<-getacc(id);if(!is.na(a)&&!is.na(accmap[a]))return(unname(accmap[a]));ck<-nk(gsub("GCST[0-9]+|EFO_?[0-9]+|\\.h\\.tsv\\.gz|\\.h\\.|_?k1$|[0-9]+$","",id));if(!is.na(nmap[ck]))return(unname(nmap[ck]));for(o in names(ov))if(grepl(nk(o),nk(id),fixed=T))return(ov[[o]]);id}
abbr<-c("Coronary Artery Disease"="CAD","Benign Prostatic Hyperplasia"="BPH","Borderline Personality Disorder"="BPD","Obstructive Sleep Apnea"="OSA","Gestational Diabetes"="GDM","Heart Failure"="HF","Myocardial Infarction"="MI","Abnormal Thrombosis"="Thrombosis","Nephrotic Syndrome"="NS","Parkinson's Disease"="PD","Major Depressive Disorder"="MDD","Inflammatory Bowel Disease"="IBD","Atrial Fibrillation"="AF","Systemic Lupus Erythematosus"="SLE","Type 1 Diabetes"="T1D","Alzheimer's Disease"="AD","Schizophrenia"="SCZ","High Myopia (retinal detachment)"="High myopia","Hay Fever"="Hay fever","Metabolic Syndrome"="MetS","Cardiometabolic Multimorbidity"="CMM","Osteoarthritis"="OA","Idiopathic Pulmonary Fibrosis"="IPF","Hearing Loss"="HL","Kidney Stones"="KS","Blood Clot"="BC","Raynaud's Disease"="Raynaud","Celiac Disease"="CeD","Crohn's Disease"="CD","Hypertrophic Cardiomyopathy"="HCM","Open-Angle Glaucoma"="Glaucoma","Chronic Obstructive Pulmonary Disease"="COPD","Lewy Body Dementia"="LBD","Pulmonary Embolism"="PE","Osteoporosis"="Osteoporosis","Restless Legs Syndrome"="RLS","Restless Leg Syndrome"="RLS")
tc<-function(s)vapply(s,function(z){w<-strsplit(z," ")[[1]];paste(vapply(w,function(x)if(nchar(x)>0&&x==toupper(x))x else paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x))),character(1)),collapse=" ")},character(1))
disp<-function(v){x<-sapply(v,mapname);tc(ifelse(x%in%names(abbr),abbr[x],x))}
co2<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[,g:=strip(gene)];co2<-co2[chr%in%as.character(1:22)];chrm<-setNames(as.integer(co2$chr),co2$g);posm<-setNames(as.numeric(co2$pos),co2$g)
POSC<-"#E08214";NEGC<-"#8073AC"
pipe<-function(file,sex){d<-fread(file,select=c("gene","disease","b1_at_maxZ","padj_acat_1"))[,g:=strip(gene)][is.finite(b1_at_maxZ)]
 d<-d[!(!is.na(chrm[g])&chrm[g]==6&posm[g]>25e6&posm[g]<34e6)]
 DICT<-fread(DICTP); d<-d[disease%in%DICT$acat_label]; d[,ds:=sign(b1_at_maxZ)]
 rhoD<-setNames(as.numeric(if(sex=="female")DICT$rhoGE_female else DICT$rhoGE_male),DICT$acat_label); d[,rho:=rhoD[disease]]
 dz<-unique(d[is.finite(rho)]$disease)
 rmap<-setNames(unique(d[disease%in%dz,.(disease,rho)])$rho,unique(d[disease%in%dz,.(disease,rho)])$disease)
 sig<-d[padj_acat_1<0.05&disease%in%dz,.(g,disease,ds)];SMx<-dcast(sig,g~disease,value.var="ds");Mx<-as.matrix(SMx[,-1]);rownames(Mx)<-SMx$g;sc<-colnames(Mx)
 out<-list();for(i in 1:(length(dz)-1))for(j in (i+1):length(dz)){a<-dz[i];b<-dz[j];if(!(a%in%sc&&b%in%sc))next;v<-Mx[,c(a,b)];v<-v[complete.cases(v),,drop=F];if(nrow(v)<15)next
   out[[length(out)+1]]<-data.table(a,b,ns=nrow(v),revrate=mean(v[,1]!=v[,2]))}
 list(P=rbindlist(out),rmap=rmap)}
saddle_tri<-function(P,rmap){
 present<-unique(rbind(P[,.(a,b)],P[,.(a=b,b=a)]));fr<-function(D)nrow(present[a%in%D&b%in%D])/(length(D)*(length(D)-1))
 gd<-rmap[is.finite(rmap)];alld<-intersect(unique(c(P$a,P$b)),names(gd))
 posd<-names(sort(gd[alld][gd[alld]>0],decreasing=T));negd<-names(sort(gd[alld][gd[alld]<0]));sel<-c(head(negd,8),head(posd,14))
 D<-sel;repeat{if(fr(D)>=.88||length(D)<15)break;pr<-present[a%in%D&b%in%D];deg<-sapply(D,function(x)sum(pr$a==x));D<-setdiff(D,names(which.min(setNames(deg,D))))}
 ordX<-D[order(gd[D])];n<-length(ordX);nneg<-sum(gd[ordX]<0);npos<-n-nneg;LS<-max(1.3,1.9*min(1,21/n))
 xcol<-setNames(n:1,ordX);yrow<-setNames(1:n,ordX);rng<-range(gd[ordX]);maxabs<-max(abs(rng));SX<-2.4/maxabs;SY<-2.4/maxabs;GAP<-.8
 xb0<-.5-GAP-maxabs*SX;yb0<-.5-GAP-maxabs*SY
 rev2<-rbind(P[,.(i=a,j=b,revrate,ns)],P[,.(i=b,j=a,revrate,ns)])[i%in%ordX&j%in%ordX]
 TT<-as.data.table(expand.grid(i=ordX,j=ordX,stringsAsFactors=F))[i!=j];TT<-merge(TT,rev2,by=c("i","j"),all.x=T);TT[,`:=`(x=xcol[i],y=yrow[j],enr=revrate-.5)];TT<-TT[x+y<=n]
 FL<-.5;FB<-signif(FL*.8,1);U<-P[a%in%ordX&b%in%ordX];U[,opp:=sign(gd[a])!=sign(gd[b])];medS<-100*median(U[opp==F]$revrate);medO<-100*median(U[opp==T]$revrate)
 LB<-data.table(d=ordX,y=yrow[ordX],v=gd[ordX]);LB[,`:=`(x2=xb0-v*SX,col=fifelse(v>0,POSC,NEGC))];vx<-npos+.5;hy<-nneg+.5;lx<-xb0-maxabs*SX-.35;by<- -0.1
 rn<-data.table(y=yrow[ordX],lab=disp(ordX));cn<-data.table(x=xcol[ordX],lab=disp(ordX))
 axy<-n+.55;cand<-c(.02,.05,.1,.15,.2,.25,.3,.4);X<-suppressWarnings(max(cand[cand<=maxabs]));if(!is.finite(X))X<-signif(maxabs*.9,1);AX<-data.table(t=c(-X,0,X));AX[,x:=xb0-t*SX]
 ML<-max(nchar(disp(ordX)));LP<-max(3.6,0.30*ML*(LS/1.9));BP<-max(2.6,0.26*ML*(LS/1.9))
 ggplot()+geom_tile(data=TT,aes(x,y,fill=enr,alpha=ns),color="white",linewidth=.2)+
  scale_fill_gradient2(low="#2166AC",mid="grey96",high="#B2182B",midpoint=0,limits=c(-FL,FL),oob=scales::squish,na.value="grey88",name="opposite-effect\nvs 50%",breaks=c(-FB,0,FB),labels=c(sprintf("-%.2g",FB),"0",sprintf("+%.2g",FB)))+
  scale_alpha_continuous(range=c(.45,1),guide="none",na.value=1)+
  annotate("segment",x=vx,xend=vx,y=.5,yend=n-npos+.5,linewidth=.45,color="grey15")+annotate("segment",y=hy,yend=hy,x=.5,xend=n-nneg+.5,linewidth=.45,color="grey15")+
  geom_rect(data=LB,aes(xmin=pmin(xb0,x2),xmax=pmax(xb0,x2),ymin=y-.42,ymax=y+.42),fill=LB$col)+annotate("segment",x=xb0,xend=xb0,y=.5,yend=n+.5,linewidth=.3,linetype=2,color="grey55")+
  annotate("segment",x=xb0-X*SX,xend=xb0+X*SX,y=axy,yend=axy,linewidth=.35,color="grey30")+geom_segment(data=AX,aes(x=x,xend=x,y=axy,yend=axy+.34),linewidth=.35,color="grey30")+
  geom_text(data=AX,aes(x=x,y=axy+1.15,label=ifelse(t==0,"0",formatC(t,format="f",digits=2))),size=1.6,color="grey30")+annotate("text",x=xb0,y=n+2.9,label="rho[GE]",parse=TRUE,size=2.1)+
  geom_text(data=rn,aes(x=lx,y=y,label=lab),hjust=1,size=LS)+geom_text(data=cn,aes(x=x,y=by,label=lab),angle=45,hjust=1,vjust=1,size=LS)+
  annotate("point",x=n*.30,y=n+2.6,shape=15,color="#C0392B",size=1.9)+annotate("text",x=n*.30+.5,y=n+2.6,parse=TRUE,label="'fertility-discordant disease pairs  ('*rho[i]%.%rho[j]*' < 0)'",hjust=0,size=1.7,color="#7B241C")+
  annotate("text",x=n*.30+.5,y=n+1.75,parse=TRUE,label=sprintf("'%.0f%% opposite gene effect'",medO),hjust=0,size=1.7,color="#7B241C")+
  annotate("point",x=n*.30,y=n+.55,shape=15,color="#2C6FBB",size=1.9)+annotate("text",x=n*.30+.5,y=n+.55,parse=TRUE,label="'fertility-concordant disease pairs  ('*rho[i]%.%rho[j]*' > 0)'",hjust=0,size=1.7,color="#1A5276")+
  annotate("text",x=n*.30+.5,y=n-.3,parse=TRUE,label=sprintf("'%.0f%% opposite gene effect'",medS),hjust=0,size=1.7,color="#1A5276")+
  annotate("point",x=n*.30,y=n-1.5,shape=15,color=POSC,size=1.9)+annotate("text",x=n*.30+.5,y=n-1.5,parse=TRUE,label="rho[GE]>0*':  positive disease-fertility'",hjust=0,size=1.7,color="grey30")+
  annotate("text",x=n*.30+.5,y=n-2.35,label="expression-effect correlation",hjust=0,size=1.7,color="grey30")+annotate("point",x=n*.30,y=n-3.55,shape=15,color=NEGC,size=1.9)+
  annotate("text",x=n*.30+.5,y=n-3.55,parse=TRUE,label="rho[GE]<0*':  negative disease-fertility'",hjust=0,size=1.7,color="grey30")+annotate("text",x=n*.30+.5,y=n-4.4,label="expression-effect correlation",hjust=0,size=1.7,color="grey30")+
  coord_fixed(ratio=1,xlim=c(lx-LP,n+4),ylim=c(by-BP,n+3.3),clip="off")+
  theme_void(base_size=7)+theme(legend.position="inside",legend.position.inside=c(.75,.60),legend.justification=c(0,1),legend.box="vertical",legend.key.size=unit(2.6,"mm"),legend.title=element_text(size=5.3),legend.text=element_text(size=4.8),legend.spacing.y=unit(.4,"mm"),plot.margin=margin(3,0,3,4))}
binGlobal<-function(P,rmap,nb){P[,gg:=rmap[a]*rmap[b]];P<-P[is.finite(gg)&gg!=0]
 P[,grp:=as.integer(cut(gg,quantile(gg,seq(0,1,length.out=nb+1)),include.lowest=T))];P[,grpf:=factor(grp,levels=1:nb)]
 bs<-P[,.(mg=mean(gg)),by=grp];P[,disc:=fifelse(bs$mg<0,"disc","conc")[match(grp,bs$grp)]];if(nb%%2==1)P[grp==(nb+1)/2,disc:="neutral"]
 H<-P[,.(med=median(revrate),n=.N),by=grp][order(grp)];H[,fit:=predict(lm(med~grp))];lb<-rep("''",nb);lb[1]<-"atop('discordant','(strong)')";lb[nb]<-"atop('concordant','(strong)')";if(nb%%2==1)lb[(nb+1)/2]<-"rho[i]*rho[j]%~~%0"
 ggplot(P,aes(grpf,revrate))+geom_hline(yintercept=.5,linetype="dotted",color="grey65")+geom_violin(aes(fill=disc),color=NA,alpha=.22,scale="width",width=.85)+
  geom_jitter(aes(size=ns,color=disc),width=.17,alpha=.14,stroke=0)+geom_boxplot(width=.14,outlier.shape=NA,linewidth=.34,fill="white",color="grey30")+geom_line(data=H,aes(grp,fit),color="grey20",linewidth=.7,inherit.aes=FALSE)+geom_point(data=H,aes(grp,med),shape=23,size=1.8,fill="white",color="black",stroke=.6)+
  scale_fill_manual(values=c(disc="#C0392B",conc="#2C6FBB",neutral="grey60"),guide="none")+scale_color_manual(values=c(disc="#C0392B",conc="#2C6FBB",neutral="grey55"),guide="none")+scale_size_continuous(range=c(.25,1.7),name="shared\ngenes")+
  scale_y_continuous(labels=function(x)x*100,limits=c(0,1))+scale_x_discrete(labels=parse(text=lb))+
  labs(x=expression("fertility coupling of disease pairs   ("*rho[i]%.%rho[j]*":  discordant "%->%" concordant)"),y="opposite-effect proportion (% shared genes)")+
  theme_classic(base_size=7)+theme(axis.text.x=element_text(size=5.7),axis.text.y=element_text(size=6),axis.title.y=element_text(size=6.5),axis.title.x=element_text(size=6,margin=margin(t=3)),legend.position="inside",legend.position.inside=c(.80,.96),legend.direction="horizontal",legend.title=element_text(size=5,vjust=.9),legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=4.6),legend.margin=margin(0,0,0,0),plot.margin=margin(18,16,10,0))}
## MAIN Fig4 D/E = MALE (consistent with panel A male lollipop); female -> supplement
Rm<-pipe("ACAT_male_gene_disease_all.csv","male");pd<-saddle_tri(Rm$P,Rm$rmap);pe<-binGlobal(copy(Rm$P),Rm$rmap,5)
Rf<-pipe("ACAT_female_gene_disease_all.csv","female");pdF<-saddle_tri(Rf$P,Rf$rmap);peF<-binGlobal(copy(Rf$P),Rf$rmap,5)

## ================= NEW PANEL: prevalence scatter (Fig 2C), NO band =================
PT<-"${ROOT}/rg_prevalence_analysis_2026-07-16/tables"
rhoP<-fread(file.path(PT,"rhoGE_finngen_mhcx.tsv"))
prevP<-fread(file.path(PT,"finngen_rg_vs_finregistry_prevalence.tsv"))[,.(phenocode,lp=log10_finregistry_prev_pct)]
dP<-merge(rhoP,prevP,by="phenocode")
mP<-melt(dP[,.(phenocode,lp,Male=rho_father,Female=rho_mother)],id.vars=c("phenocode","lp"),variable.name="axis",value.name="rho")[is.finite(rho)&is.finite(lp)]
mP[,axis:=factor(axis,levels=c("Male","Female"))]
sP<-mP[,{ct<-cor.test(rho,lp);.(lab=sprintf("'%s: r=%.2f, '*italic(P)==%s*' (n=%d)'",axis,ct$estimate,pmP(ct$p.value),.N))},by=axis]
mP[,prevpct:=10^lp]
pPrev<-ggplot(mP,aes(rho,prevpct,color=axis))+geom_vline(xintercept=0,color="grey85",linetype=2,linewidth=.3)+
  geom_point(size=.5,alpha=.28)+geom_smooth(aes(linetype=axis),method="lm",formula=y~x,se=FALSE,linewidth=.8)+
  geom_text(data=sP[axis=="Male"],aes(x=-Inf,y=Inf,label=lab),color=COLS["Male"],hjust=-.05,vjust=1.5,size=2.3,parse=TRUE,inherit.aes=FALSE)+
  geom_text(data=sP[axis=="Female"],aes(x=-Inf,y=Inf,label=lab),color=COLS["Female"],hjust=-.05,vjust=3,size=2.3,parse=TRUE,inherit.aes=FALSE)+
  scale_color_manual(values=COLS,name=NULL)+scale_linetype_manual(values=LTS,name=NULL)+guides(linetype="none")+
  scale_y_log10(breaks=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))+
  labs(x=expression("disease-fertility "*rho[GE]*" (expression effect)"),y="Population prevalence (%)")+th+
  theme(legend.position=c(.99,.02),legend.justification=c(1,0),legend.background=element_rect(fill=alpha("white",.7),color=NA))

## ---- prevalence violin (old Fig2B) -> SUPPLEMENT ----
mP[,grp:=factor(fifelse(rho>0,"rho>0","rho<0"),levels=c("rho<0","rho>0"))]
stV<-mP[,{w<-wilcox.test(lp~grp);.(p=w$p.value,y=max(lp)+.3)},by=axis]; cat("prevalence violin Mann-Whitney:\n");print(stV)
set.seed(1)
medP<-mP[,.(m=median(prevpct)),by=.(axis,grp)]
pPrevB<-ggplot(mP,aes(grp,prevpct))+geom_violin(aes(fill=grp),color=NA,alpha=.24,scale="width",width=.85)+
  geom_jitter(aes(color=grp),width=.28,height=0,size=.28,alpha=.16,stroke=0)+
  geom_boxplot(width=.16,outlier.shape=NA,linewidth=.35,fill="white",color="grey30")+
  geom_point(data=medP,aes(grp,m),shape=23,size=1.5,fill="white",color="black",stroke=.5)+facet_wrap(~axis)+
  scale_fill_manual(values=c(`rho<0`=unname(COLS["Male"]),`rho>0`=unname(COLS["Female"])),guide="none")+
  scale_color_manual(values=c(`rho<0`=unname(COLS["Male"]),`rho>0`=unname(COLS["Female"])),guide="none")+
  scale_x_discrete(labels=c(`rho<0`=expression(rho[GE]<0),`rho>0`=expression(rho[GE]>0)))+
  geom_text(data=stV,aes(x=1.5,y=10^y,label=sprintf("Mann-Whitney~italic(P)==%s",pmP(p))),parse=TRUE,size=2.3,inherit.aes=FALSE)+
  scale_y_log10(breaks=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))+
  labs(x=expression("disease-fertility "*rho[GE]),y="Population prevalence (%)")+
  theme_bw(BS)+theme(panel.grid=element_blank(),strip.background=element_rect(fill="grey92",color=NA),strip.text=element_text(face="bold",size=BS-1),axis.text=element_text(size=BS-1,color="grey20"))
ggsave(file.path(MY,"fusio2trait/figures_final/Fig_prevalence_violin_SUPP.pdf"),pPrevB,width=130,height=85,units="mm",device=cairo_pdf,bg="white")

## ================= NEW PANEL: 4-class distribution (Fig 4A), male+female =================
d4<-CF[,.N,by=.(Sex=fifelse(fitness=="male","Male","Female"),cc=combo)][,prop:=N/sum(N),by=Sex]
lev4<-c("++","--","+-","-+"); lab4<-c("++"="d+f+","--"="d-f-","+-"="d+f-","-+"="d-f+")
cols4<-c("++"="#C0392B","--"="#E67E22","+-"="#2C5985","-+"="#5DADE2")
d4[,cc:=factor(cc,levels=lev4)];d4[,Sex:=factor(Sex,levels=c("Male","Female"))]
p4c<-ggplot(d4,aes(cc,prop,fill=cc))+geom_col(width=.72)+geom_text(aes(label=N),vjust=-.3,size=1.9)+
  facet_wrap(~Sex)+scale_fill_manual(values=cols4,labels=lab4,name=NULL)+scale_x_discrete(labels=lab4)+
  scale_y_continuous(expand=expansion(mult=c(0,.14)),labels=scales::percent)+
  labs(x="effect-direction class",y="proportion of gene-disease pairs")+th+
  theme(legend.position="none",axis.text.x=element_text(size=BS-2),strip.background=element_rect(fill="grey92",color=NA),strip.text=element_text(face="bold",size=BS-1))
cat("4-class proportions by sex:\n");print(dcast(d4,Sex~cc,value.var="prop")[,lapply(.SD,function(x)if(is.numeric(x))round(100*x,1) else x)])

## ================= ASSEMBLE Fig 2 (A lollipop / B concordance / C prevalence) =================
## standalone panel saves (2026-08-14: Fig2 rebuilt as A + B|C in LaTeX)
ggsave(file.path(MY,"fusio2trait/figures_final/Fig2_panelA.pdf"),pa,width=183,height=78,units="mm",device=cairo_pdf,bg="white")
ggsave(file.path(MY,"fusio2trait/figures_final/Fig2_panelB_conc.pdf"),pb,width=89,height=80,units="mm",device=cairo_pdf,bg="white")
ggsave(file.path(MY,"fusio2trait/figures_final/Fig2_panelPrev.pdf"),pPrev,width=89,height=80,units="mm",device=cairo_pdf,bg="white")
row2b<-plot_grid(pb,pPrev,ncol=2,rel_widths=c(1,1),labels=c("B","C"),label_size=12,label_fontface="bold")
fig2<-plot_grid(pa,row2b,ncol=1,rel_heights=c(0.92,1.0),labels=c("A",""),label_size=12,label_fontface="bold")
ggsave(file.path(MY,"Fig2_3panel.pdf"),fig2,width=183,height=172,units="mm",device=cairo_pdf,bg="white")
ggsave(file.path(MY,"Fig2_3panel.png"),fig2,width=183,height=172,units="mm",dpi=190,bg="white")

## ================= ASSEMBLE Fig 4 (2x2: A 4class / B %antag / C triangle / D violin) =================
r4a<-plot_grid(p4c,pc,ncol=2,rel_widths=c(1,1),labels=c("A","B"),label_size=12,label_fontface="bold")
## export the individual panels so the figures can be re-assembled by layer
saveRDS(list(fourclass=p4c, antag=pc, matrix=pd, bins=pe),
        file.path(MY,"fusios_pipeline","fig4_panels.rds"))
r4b<-plot_grid(pd,pe,ncol=2,rel_widths=c(1.15,1),labels=c("C","D"),label_size=12,label_fontface="bold")
fig4<-plot_grid(r4a,r4b,ncol=1,rel_heights=c(1,1))
ggsave(file.path(MY,"Fig4_4panel.pdf"),fig4,width=183,height=158,units="mm",device=cairo_pdf,bg="white")
ggsave(file.path(MY,"Fig4_4panel.png"),fig4,width=183,height=158,units="mm",dpi=190,bg="white")

## ---- FEMALE supplement (triangle + violin) unchanged ----
rowF<-plot_grid(pdF,peF,ncol=2,rel_widths=c(1.15,1),labels=c("A","B"),label_size=12,label_fontface="bold")
ggsave(file.path(MY,"fusio2trait/figures_final/Fig_intralocus_SUPP_female.pdf"),rowF,width=183,height=95,units="mm",device=cairo_pdf,bg="white")
cat("DONE Fig2_3panel + Fig4_4panel. prevalence panel numbers:\n");print(sP)
cat("panel c weighted r/P:\n");print(rBL[,.(sex,r,p=signif(p,2))])
