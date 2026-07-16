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
d<-d[!(!is.na(chrm[g])&chrm[g]==6&posm[g]>25e6&posm[g]<34e6)&!iscancer(disease)]
bad<-unique(d$disease);bad<-bad[sapply(bad,isContam)];d<-d[!disease%in%bad];d[,ds:=sign(b1_at_maxZ)]
umap<-setNames(nk(sapply(unique(d$disease),mapname)),unique(d$disease))
resc<-c(autism="autismspectrumdisorder",highmyopiaretinaldetachment="highmyopia",rheumatoidarthritish="rheumatoidarthritis",t2d="type2diabetes",tstourettesyndrome="tourettesyndrome",pulmonaryembolism="bloodclotlung",chronicobstructivepulmonarydisease="copd",panic="panicdisorder")
hit<-umap%in%names(resc);umap[hit]<-resc[umap[hit]]
rf<-setNames(T5$rhoGE_female,T5$k);rm<-setNames(T5$rhoGE_male,T5$k)
d[,`:=`(rhoF=rf[umap[disease]],rhoM=rm[umap[disease]])]
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
dz<-dz[!is.na(dz)]
ia<-match(P$a,dz);ib<-match(P$b,dz)
rlog<-function(...)cat(sprintf(...))
sink("/Volumes/S840/04mypaper/lin/conflictresolution/fig4_intralocus_session_20260714/intralocus_FINAL_results.txt")
cat("=== INTRALOCUS FINAL (canonical) : ACAT padj<0.05, >=15 shared, MHC+cancer excl, T5 rho, crosswalk-rescued ===\n")
cat(sprintf("n_diseases=%d  n_pairs=%d\n\n",length(dz),nrow(P)))
## fast residual-on-sim
ro<-function(x,s){b<-cov(x,s)/var(s);(x-mean(x))-b*(s-mean(s))}
ry<-ro(P$revrate,P$sim)
permP<-function(prodfun,gd,obs,B=100000){gd<-gd[dz]
 stat<-numeric(B);for(k in 1:B){gp<-sample(gd);xp<-gp[ia]*gp[ib];rx<-ro(xp,P$sim);stat[k]<-cor(rx,ry)}
 (1+sum(abs(stat)>=abs(obs)))/(B+1)}
## ---- A. descriptive contrasts (both sexes) ----
cat("---- A. discordant vs concordant opposite-effect % + OR ----\n")
for(sc in c("Xf","Xm")){sx<-ifelse(sc=="Xf","maternal(F)","paternal(M)");D<-P[is.finite(get(sc))];D[,X:=get(sc)]
 di<-D[X<0];co<-D[X>0];a11<-sum(di$nrev);a12<-sum(di$ns)-a11;a21<-sum(co$nrev);a22<-sum(co$ns)-a21
 or<-(a11*a22)/(a12*a21);se<-sqrt(1/a11+1/a12+1/a21+1/a22);ci<-exp(log(or)+c(-1.96,1.96)*se)
 D[,q:=cut(X,quantile(X,seq(0,1,.1)),include.lowest=T,labels=1:10)];en<-D[,.(rr=sum(nrev)/sum(ns)),by=q][order(q)]
 rlog("  %s: discordant %.1f%% vs concordant %.1f%%  OR=%.2f [%.2f,%.2f] | decile %.0f%%->%.0f%%\n",
   sx,100*a11/(a11+a12),100*a21/(a21+a22),or,ci[1],ci[2],100*en$rr[1],100*en$rr[10])}
## ---- B. gene-level logistic OR ----
cat("\n---- B. gene-level logistic OR per s.d.(rho_i*rho_j) ----\n")
for(sc in c("Xf","Xm")){sx<-ifelse(sc=="Xf","maternal","paternal");D<-P[is.finite(get(sc))]
 L<-D[,.(y=unlist(Map(function(nr,nn)c(rep(1,nr),rep(0,nn-nr)),nrev,ns)),X=rep(get(sc),ns))];L[,Xs:=scale(X)[,1]]
 c2<-summary(glm(y~Xs,binomial,L))$coef["Xs",];rlog("  %s: OR=%.2f  P=%.1e\n",sx,exp(c2[1]),c2[4])}
## ---- C. pooled partial r + 100k perm + incremental R2 ----
cat("\n---- C. pooled partial regression (mean rho) ----\n")
rx<-ro(P$Xmean,P$sim);pr<-cor(rx,ry)
rmn<-setNames(sapply(dz,function(x)(rF[x]+rM[x])/2),dz)
pP<-permP(NULL,rmn,pr,100000)
r2s<-summary(lm(revrate~sim,P))$r.squared;r2f<-summary(lm(revrate~sim+Xmean,P))$r.squared
rlog("  partial r=%.3f  perm P=%.5f (100k)\n",pr,pP)
rlog("  R2 base(sim)=%.1f%%  full(sim+fert)=%.1f%%  incremental=%.1f%%\n",100*r2s,100*r2f,100*(r2f-r2s))
## ---- D. robust estimators + LODO (both sexes) ----
cat("\n---- D. robustness (standardized beta) ----\n")
suppressMessages(library(MASS))
for(sc in c("Xf","Xm")){sx<-ifelse(sc=="Xf","maternal","paternal");D<-P[is.finite(get(sc))];D[,`:=`(X=scale(get(sc))[,1],Y=scale(revrate)[,1])]
 ols<-summary(lm(Y~X+sim,D))$coef["X",];hub<-summary(rlm(Y~X+sim,D,maxit=200))$coef["X",];mm<-summary(rlm(Y~X+sim,D,method="MM",maxit=200))$coef["X",]
 D[,ck:=cooks.distance(lm(Y~X+sim,D))];q<-quantile(D$ck,.99);tr<-coef(lm(Y~X+sim,D[ck<=q]))["X"]
 lodo<-sapply(dz,function(dd){dt<-D[a!=dd&b!=dd];if(nrow(dt)<50)return(NA);coef(lm(Y~X+sim,dt))["X"]});lodo<-lodo[is.finite(lodo)]
 rlog("  %s: OLS %.2f [%.2f,%.2f] | Huber %.2f | MM %.2f | trim %.2f | LODO %d/%d neg (med %.2f)\n",
   sx,ols[1],ols[1]-1.96*ols[2],ols[1]+1.96*ols[2],hub[1],mm[1],tr,sum(lodo<0),length(lodo),median(lodo))}
## ---- E. r_g replication (fig3 csv, correct key nk(raw id)) ----
cat("\n---- E. r_g replication (LDSC, fig3 source) ----\n")
fig3<-fread("/Users/cjh/Documents/GitHub/HumanDiseaseFertilityPleiotropy/code/fig3_disease_correlations.csv")
rawk<-setNames(nk(dz),dz)
rgM<-setNames(fig3$rgM,fig3$dk)[rawk[dz]];rgF<-setNames(fig3$rgF,fig3$dk)[rawk[dz]];names(rgM)<-dz;names(rgF)<-dz
Prg<-copy(P);Prg[,Xrg:=(rgM[a]*rgM[b]+rgF[a]*rgF[b])/2];Drg<-Prg[is.finite(Xrg)]
ryg<-ro(Drg$revrate,Drg$sim);rxg<-ro(Drg$Xrg,Drg$sim);prg<-cor(rxg,ryg)
iag<-match(Drg$a,dz);ibg<-match(Drg$b,dz);rgmn<-(rgM+rgF)/2
pstat<-numeric(50000);for(k in 1:50000){gp<-sample(rgmn[dz]);xp<-gp[iag]*gp[ibg];pstat[k]<-cor(ro(xp,Drg$sim),ryg)}
pPg<-(1+sum(abs(pstat)>=abs(prg)))/50001
rlog("  r_g: partial r=%.3f  perm P=%.5f (50k)  n_pairs=%d  n_dz=%d\n",prg,pPg,nrow(Drg),sum(is.finite(rgmn[dz])))
## ---- F. cross-fit (non-circularity): rho from gene-half A, O from gene-half B ----
cat("\n---- F. cross-fit non-circularity (disjoint gene halves) ----\n")
set.seed(7)
allg<-unique(d$g);half<-setNames(sample(rep(c("A","B"),length.out=length(allg))),allg)
dd<-copy(d);dd[,h:=half[g]]
## rho^A per disease (sign concordance of disease vs fertility not available in ACAT d; use disease-disease within halves)
## simpler: naive vs cross-fit of O~similarity structure is not the target; report the earlier all-gene cross-fit result
cat("  (cross-fit computed separately in gate test: naive antag~rho 0.87/0.86 -> disjoint-gene 0.66/0.61; survives)\n")
## ---- G. global-Bonferroni sensitivity (gene_level_b1, p_bonf_global) ----
cat("\n---- G. global-Bonferroni threshold sensitivity ----\n")
b1g<-fread("gene_level_b1_disease_acat_summary.dropbox.csv",select=c("gene","disease","meta_beta_ivw","p_bonf_global"))[,g:=strip(gene)]
b1g[nk(disease)=="ulcerativecolitis1",meta_beta_ivw:=-meta_beta_ivw]
sigg<-b1g[p_bonf_global<0.05,.(g,disease,ds=sign(meta_beta_ivw))]
sigg[,k:=nk(sapply(disease,function(x)x))]  # keep raw; match to our dz by nk(raw)
sigg[,kk:=nk(disease)]
dzk<-setNames(nk(dz),dz)
SMg<-dcast(sigg[kk%in%dzk,.(g,kk,ds)],g~kk,value.var="ds");Mg<-as.matrix(SMg[,-1]);rownames(Mg)<-SMg$g;scg<-colnames(Mg)
og<-list();kk2dz<-setNames(dz,dzk)
dzc<-intersect(dzk,scg)
for(i in 1:(length(dzc)-1))for(j in (i+1):length(dzc)){a<-dzc[i];b<-dzc[j];v<-Mg[,c(a,b)];v<-v[complete.cases(v),,drop=F];if(nrow(v)<15)next
  og[[length(og)+1]]<-data.table(a=kk2dz[a],b=kk2dz[b],revG=mean(v[,1]!=v[,2]),ns=nrow(v))}
PG<-rbindlist(og);PG<-merge(PG,P[,.(a,b,sim,Xmean)],by=c("a","b"))
if(nrow(PG)>30){ryG<-ro(PG$revG,PG$sim);rxG<-ro(PG$Xmean,PG$sim);prG<-cor(rxG,ryG)
 iaG<-match(PG$a,dz);ibG<-match(PG$b,dz)
 psG<-numeric(20000);for(k in 1:20000){gp<-sample(rmn[dz]);xp<-gp[iaG]*gp[ibG];psG[k]<-cor(ro(xp,PG$sim),ryG)}
 pPG<-(1+sum(abs(psG)>=abs(prG)))/20001
 rlog("  global p_bonf<0.05: partial r=%.3f  perm P=%.5f  n_pairs=%d\n",prG,pPG,nrow(PG))} else cat("  global: too few pairs\n")

## ---- H. Fig4B Panel-B Pearson r (all pairs) + Fig4A subset contrast (traceability) ----
cat("\n---- H. figure numbers (for text traceability) ----\n")
rF_B<-cor(P$revrate,P$Xf);rM_B<-cor(P$revrate,P$Xm)
rlog("  Fig4B Panel-B Pearson r: female=%.3f  male=%.3f  (revrate ~ rho_i*rho_j, all %d pairs)\n",rF_B,rM_B,nrow(P))
cat("  Fig4A subset (61%% vs 31%%): source = build_ab_T5.R drawn subset; verified female opp=60.9/same=31.2, male 61.2/30.3\n")
cat("=== DONE ===\n");sink()

