## Coordinated A+B, ONE pipeline: sex-stratified ACAT (disease z + fitness z) -> rho_GE + reversal.
## A = LOWER-TRIANGLE saddle (curated diseases; symmetric heatmap+bars -> keep left-bottom half, no info lost).
## B = 6-bin reversal gradient. Same data/metric/style. Female=main, Male=supplement.
suppressMessages({library(data.table);library(ggplot2);library(cowplot);library(readxl)});options(width=200);set.seed(1);nk<-function(x)gsub("[^a-z0-9]","",tolower(x));strip<-function(x)sub("[.].*","",x)
BD<-"/Volumes/S840/04mypaper/lin/conflictresolution/fig4_intralocus_session_20260714/";setwd("/Volumes/X10Pro/data/synergistic/mr_res")
st<-as.data.table(read_excel("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/Supplementary Table_v3.xlsx",sheet="Supplementary Table 1",skip=3));setnames(st,1:7,c("No","Trait","Type","Category","Source","Accession","PMID"));st<-st[!is.na(Trait)&grepl("^[0-9]",No)]
accmap<-setNames(st$Trait,gsub("[^A-Za-z0-9]","",st$Accession));nmap<-setNames(st$Trait,nk(st$Trait))
getacc<-function(x){m<-regmatches(x,regexpr("GCST[0-9]+",x));if(length(m))m else NA}
ov<-c(cad="Coronary Artery Disease",ibd="Inflammatory Bowel Disease",bip="Bipolar Disorder",adhd="ADHD",scz="Schizophrenia",ptsd="PTSD",Crohns="Crohn's Disease",Raynaud="Raynaud's Disease",Hyperplasiaofprostate="Benign Prostatic Hyperplasia",Highmyopia="High Myopia",osteoarthritis="Osteoarthritis",Heartfailure="Heart Failure",myocardialinfarction="Myocardial Infarction",majordepressivedisorder="Major Depressive Disorder",parkinson="Parkinson's Disease",alzheimer="Alzheimer's Disease",sle="Systemic Lupus Erythematosus",t1d="Type 1 Diabetes",atrialfibrillation="Atrial Fibrillation",endometriosis="Endometriosis",osteoporosis="Osteoporosis",hayfever="Hay Fever",bloodclot="Blood Clot",cardiometabolicmulti="Cardiometabolic Multimorbidity",hypertrophiccardiomyopathy="Hypertrophic Cardiomyopathy",metabolicsyndrome="Metabolic Syndrome",celiac="Celiac Disease",kidneystones="Kidney Stones")
idov<-c(Blood_Clot_Lung="Pulmonary Embolism",Chronic_obstructive_pulmonary_disease="Chronic Obstructive Pulmonary Disease")
mapname<-function(id){if(id%in%names(idov))return(unname(idov[id]));a<-getacc(id);if(!is.na(a)&&!is.na(accmap[a]))return(unname(accmap[a]));ck<-nk(gsub("GCST[0-9]+|EFO_?[0-9]+|\\.h\\.tsv\\.gz|\\.h\\.|_?k1$|[0-9]+$","",id));if(!is.na(nmap[ck]))return(unname(nmap[ck]));for(o in names(ov))if(grepl(nk(o),nk(id),fixed=T))return(ov[[o]]);id}
abbr<-c("Coronary Artery Disease"="CAD","Benign Prostatic Hyperplasia"="BPH","Borderline Personality Disorder"="BPD","Obstructive Sleep Apnea"="OSA","Gestational Diabetes"="GDM","Heart Failure"="HF","Myocardial Infarction"="MI","Abnormal Thrombosis"="Thrombosis","Nephrotic Syndrome"="NS","Parkinson's Disease"="PD","Major Depressive Disorder"="MDD","Inflammatory Bowel Disease"="IBD","Atrial Fibrillation"="AF","Systemic Lupus Erythematosus"="SLE","Type 1 Diabetes"="T1D","Alzheimer's Disease"="AD","Schizophrenia"="SCZ","High Myopia (retinal detachment)"="High myopia","Hay Fever"="Hay fever","Metabolic Syndrome"="MetS","Cardiometabolic Multimorbidity"="CMM","Osteoarthritis"="OA","Idiopathic Pulmonary Fibrosis"="IPF","Hearing Loss"="HL","Kidney Stones"="KS","Blood Clot"="BC","Raynaud's Disease"="Raynaud","Celiac Disease"="CeD","Crohn's Disease"="CD","Hypertrophic Cardiomyopathy"="HCM","Open-Angle Glaucoma"="Glaucoma","Chronic Obstructive Pulmonary Disease"="COPD","Lewy Body Dementia"="LBD","Pulmonary Embolism"="PE","Osteoporosis"="Osteoporosis","Restless Legs Syndrome"="RLS","Restless Leg Syndrome"="RLS")
tc<-function(s)vapply(s,function(z){w<-strsplit(z," ")[[1]];paste(vapply(w,function(x)if(nchar(x)>0&&x==toupper(x))x else paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x))),character(1)),collapse=" ")},character(1))
disp<-function(v){x<-sapply(v,mapname);tc(ifelse(x%in%names(abbr),abbr[x],x))}
iscancer<-function(x)grepl("cancer|carcinoma|lymphoma|leukemia|myeloma|melanoma|Fibroid",x,ignore.case=T)
isContam<-function(x){z<-nk(x);z=="bloodclot"||grepl("celiac|sclerosis|ankylos|covid|narcoleps|metabolicsyndrome|cardiometabolic|manic|mentalhealth",z)}   # verified non-97 (audited against Supp)
co<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[,g:=strip(gene)];co<-co[chr%in%as.character(1:22)];chrm<-setNames(as.integer(co$chr),co$g);posm<-setNames(as.numeric(co$pos),co$g)
POSC<-"#E08214";NEGC<-"#8073AC";fmt<-function(z)ifelse(z==0,"0",formatC(z,format="f",digits=2))
## AUTHORITATIVE rho_GE from Supplementary Table 5 (make_tables.R: top_beta/top_se z, |z|<9 cap, MHC-excl, n>=30) -- do NOT recompute
T5<-as.data.table(read_excel("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/Supplementary Table_v3.xlsx",sheet="Supplementary Table 5",skip=1));T5<-T5[!is.na(Disease)];T5[,k:=nk(Disease)]
pipe<-function(file,sex){d<-fread(file,select=c("gene","disease","b1_at_maxZ","padj_acat_1"))[,g:=strip(gene)][is.finite(b1_at_maxZ)]
 d<-d[!(!is.na(chrm[g])&chrm[g]==6&posm[g]>25e6&posm[g]<34e6)]
 ## AUTHORITATIVE mapping via disease_dictionary (b1_disease_id = exact ACAT label -> canonical rho_GE); fixes mapname *cancerwu drop
 DICT<-fread("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/disease_dictionary.csv")
 d<-d[disease%in%DICT$acat_label]     # restrict to authoritative 97 (drops contamination automatically)
 d[,ds:=sign(b1_at_maxZ)]
 rhoD<-setNames(as.numeric(if(sex=="female")DICT$rhoGE_female else DICT$rhoGE_male),DICT$acat_label)
 d[,rho:=rhoD[disease]]
 dz<-unique(d[is.finite(rho)]$disease)
 rmap<-setNames(unique(d[disease%in%dz,.(disease,rho)])$rho,unique(d[disease%in%dz,.(disease,rho)])$disease)
 sig<-d[padj_acat_1<0.05&disease%in%dz,.(g,disease,ds)];SMx<-dcast(sig,g~disease,value.var="ds");Mx<-as.matrix(SMx[,-1]);rownames(Mx)<-SMx$g;sc<-colnames(Mx)
 out<-list();for(i in 1:(length(dz)-1))for(j in (i+1):length(dz)){a<-dz[i];b<-dz[j];if(!(a%in%sc&&b%in%sc))next;v<-Mx[,c(a,b)];v<-v[complete.cases(v),,drop=F];if(nrow(v)<15)next
   out[[length(out)+1]]<-data.table(a,b,ns=nrow(v),revrate=mean(v[,1]!=v[,2]))}
 cat(sex,": diseases with T5 rho =",length(dz)," drawable pairs =",length(out),"\n")
 list(P=rbindlist(out),rmap=rmap)}

## ---------- LOWER-LEFT TRIANGLE saddle (Panel A) ----------
## Disease-selection criterion: keep diseases with finite rho_GE from >=30 genes (stable) and
## >=15 shared genes per drawn pair; from those take the 8 most fertility-FAVOURED (most negative rho_GE)
## + 14 most fertility-COSTLY (most positive rho_GE) to span the axis, then greedily drop the
## least-connected disease until >=88% of possible pairs have data (keeps the heatmap dense).
## Bars = ORIGINAL diverging style (rho_GE, +/- from a zero baseline) on the two triangle legs (left + bottom);
## black cross-lines split the +/- (opposite vs same-sign) blocks.
saddle_tri<-function(P,rmap,full=FALSE,sm=0){
 present<-unique(rbind(P[,.(a,b)],P[,.(a=b,b=a)]));fr<-function(D)nrow(present[a%in%D&b%in%D])/(length(D)*(length(D)-1))
 gd<-rmap[is.finite(rmap)];alld<-intersect(unique(c(P$a,P$b)),names(gd))
 posd<-names(sort(gd[alld][gd[alld]>0],decreasing=T));negd<-names(sort(gd[alld][gd[alld]<0]));sel<-c(head(negd,8),head(posd,14))
 if(full){D<-alld}else{D<-sel;repeat{if(fr(D)>=.88||length(D)<15)break;pr<-present[a%in%D&b%in%D];deg<-sapply(D,function(x)sum(pr$a==x));D<-setdiff(D,names(which.min(setNames(deg,D))))}}
 ordX<-D[order(gd[D])];n<-length(ordX);nneg<-sum(gd[ordX]<0);npos<-n-nneg;LS<-max(1.3,1.9*min(1,21/n))
 xcol<-setNames(n:1,ordX);yrow<-setNames(1:n,ordX)
 rng<-range(gd[ordX]);maxabs<-max(abs(rng));SX<-2.4/maxabs;SY<-2.4/maxabs;GAP<-.8
 xb0<-.5-GAP-maxabs*SX;yb0<-.5-GAP-maxabs*SY                        # zero baselines just outside the left / bottom edges
 rev2<-rbind(P[,.(i=a,j=b,revrate,ns)],P[,.(i=b,j=a,revrate,ns)])[i%in%ordX&j%in%ordX]
 TT<-as.data.table(expand.grid(i=ordX,j=ordX,stringsAsFactors=F))[i!=j]
 TT<-merge(TT,rev2,by=c("i","j"),all.x=T);TT[,`:=`(x=xcol[i],y=yrow[j],enr=revrate-.5)]
 TT<-TT[x+y<=n]                                                     # keep lower-left half (symmetric across anti-diagonal)
 if(sm>0){lk<-setNames(TT[is.finite(enr)]$enr,paste(TT[is.finite(enr)]$x,TT[is.finite(enr)]$y))    # 2D neighbourhood smoothing to reveal block structure
   TT[,enr:=vapply(seq_len(.N),function(ii){v<-lk[as.vector(outer((x[ii]-sm):(x[ii]+sm),(y[ii]-sm):(y[ii]+sm),function(a,b)paste(a,b)))];m<-mean(v,na.rm=TRUE);if(is.nan(m))NA_real_ else m},numeric(1))]}
 FL<-if(sm>0)max(.08,signif(quantile(abs(TT$enr),.92,na.rm=TRUE),1))else .5;FB<-signif(FL*.8,1)
 U<-P[a%in%ordX&b%in%ordX];U[,opp:=sign(gd[a])!=sign(gd[b])];medS<-100*median(U[opp==F]$revrate);medO<-100*median(U[opp==T]$revrate)
 LB<-data.table(d=ordX,y=yrow[ordX],v=gd[ordX]);LB[,`:=`(x2=xb0-v*SX,col=fifelse(v>0,POSC,NEGC))]   # left marginal: +rho left, -rho right
 BB<-data.table(d=ordX,x=xcol[ordX],v=gd[ordX]);BB[,`:=`(y2=yb0+v*SY,col=fifelse(v>0,POSC,NEGC))]    # bottom marginal: +rho up, -rho down
 vx<-npos+.5;hy<-nneg+.5;lx<-xb0-maxabs*SX-.35;by<- -0.1   # bottom bars removed -> column labels sit just under heatmap
 rn<-data.table(y=yrow[ordX],lab=disp(ordX));cn<-data.table(x=xcol[ordX],lab=disp(ordX))
 axy<-n+.55;cand<-c(.02,.05,.1,.15,.2,.25,.3,.4);X<-suppressWarnings(max(cand[cand<=maxabs]));if(!is.finite(X))X<-signif(maxabs*.9,1);AX<-data.table(t=c(-X,0,X));AX[,x:=xb0-t*SX]   # rho_GE bar axis: 3 nice ticks
 ML<-max(nchar(disp(ordX)));LP<-max(3.6,0.30*ML*(LS/1.9));BP<-max(2.6,0.26*ML*(LS/1.9))   # adaptive room for long row(LP)/angled-column(BP) labels
 ggplot()+
  geom_tile(data=TT,aes(x,y,fill=enr,alpha=ns),color="white",linewidth=.2)+
  scale_fill_gradient2(low="#2166AC",mid="grey96",high="#B2182B",midpoint=0,limits=c(-FL,FL),oob=scales::squish,na.value="grey88",name="opposite-effect\nvs 50%",breaks=c(-FB,0,FB),labels=c(sprintf("-%.2g",FB),"0",sprintf("+%.2g",FB)))+
  scale_alpha_continuous(range=c(.45,1),guide="none",na.value=1)+
  annotate("segment",x=vx,xend=vx,y=.5,yend=n-npos+.5,linewidth=.45,color="grey15")+          # +/- column divider (opposite vs same block)
  annotate("segment",y=hy,yend=hy,x=.5,xend=n-nneg+.5,linewidth=.45,color="grey15")+           # +/- row divider
  geom_rect(data=LB,aes(xmin=pmin(xb0,x2),xmax=pmax(xb0,x2),ymin=y-.42,ymax=y+.42),fill=LB$col)+
  annotate("segment",x=xb0,xend=xb0,y=.5,yend=n+.5,linewidth=.3,linetype=2,color="grey55")+
  annotate("segment",x=xb0-X*SX,xend=xb0+X*SX,y=axy,yend=axy,linewidth=.35,color="grey30")+
  geom_segment(data=AX,aes(x=x,xend=x,y=axy,yend=axy+.34),linewidth=.35,color="grey30")+
  geom_text(data=AX,aes(x=x,y=axy+1.15,label=ifelse(t==0,"0",formatC(t,format="f",digits=2))),size=1.6,color="grey30")+
  annotate("text",x=xb0,y=n+2.9,label="rho[GE]",parse=TRUE,size=2.1)+
  geom_text(data=rn,aes(x=lx,y=y,label=lab),hjust=1,size=LS)+
  geom_text(data=cn,aes(x=x,y=by,label=lab),angle=45,hjust=1,vjust=1,size=LS)+
  annotate("point",x=n*.30,y=n+2.6,shape=15,color="#C0392B",size=1.9)+
  annotate("text",x=n*.30+.5,y=n+2.6,parse=TRUE,label="'fertility-discordant disease pairs  ('*rho[i]%.%rho[j]*' < 0)'",hjust=0,size=1.7,color="#7B241C")+
  annotate("text",x=n*.30+.5,y=n+1.75,parse=TRUE,label=sprintf("'%.0f%% opposite gene effect'",medO),hjust=0,size=1.7,color="#7B241C")+
  annotate("point",x=n*.30,y=n+.55,shape=15,color="#2C6FBB",size=1.9)+
  annotate("text",x=n*.30+.5,y=n+.55,parse=TRUE,label="'fertility-concordant disease pairs  ('*rho[i]%.%rho[j]*' > 0)'",hjust=0,size=1.7,color="#1A5276")+
  annotate("text",x=n*.30+.5,y=n-.3,parse=TRUE,label=sprintf("'%.0f%% opposite gene effect'",medS),hjust=0,size=1.7,color="#1A5276")+
  annotate("point",x=n*.30,y=n-1.5,shape=15,color=POSC,size=1.9)+
  annotate("text",x=n*.30+.5,y=n-1.5,parse=TRUE,label="rho[GE]>0*':  positive disease-fertility'",hjust=0,size=1.7,color="grey30")+
  annotate("text",x=n*.30+.5,y=n-2.35,label="expression-effect correlation",hjust=0,size=1.7,color="grey30")+
  annotate("point",x=n*.30,y=n-3.55,shape=15,color=NEGC,size=1.9)+
  annotate("text",x=n*.30+.5,y=n-3.55,parse=TRUE,label="rho[GE]<0*':  negative disease-fertility'",hjust=0,size=1.7,color="grey30")+
  annotate("text",x=n*.30+.5,y=n-4.4,label="expression-effect correlation",hjust=0,size=1.7,color="grey30")+
  coord_fixed(ratio=1,xlim=c(lx-LP,n+4),ylim=c(by-BP,n+3.3),clip="off")+
  theme_void(base_size=7)+theme(legend.position="inside",legend.position.inside=c(.75,.60),legend.justification=c(0,1),legend.box="vertical",legend.key.size=unit(2.6,"mm"),legend.title=element_text(size=5.3),legend.text=element_text(size=4.8),legend.spacing.y=unit(.4,"mm"),plot.margin=margin(3,0,3,4))}
## GLOBAL binning: one quantile cut over the whole signed coupling axis gg=rho_i*rho_j (nb bins),
## running strong-discordant(most negative) -> strong-concordant(most positive). Single monotone series,
## violin/box coloured by bin-mean sign, jittered pairs (size=shared genes), white median diamonds.
binGlobal<-function(P,rmap,nb){P[,gg:=rmap[a]*rmap[b]];P<-P[is.finite(gg)&gg!=0]
 P[,grp:=as.integer(cut(gg,quantile(gg,seq(0,1,length.out=nb+1)),include.lowest=T))]
 P[,grpf:=factor(grp,levels=1:nb)]
 bs<-P[,.(mg=mean(gg)),by=grp];P[,disc:=fifelse(bs$mg<0,"disc","conc")[match(grp,bs$grp)]]
 if(nb%%2==1)P[grp==(nb+1)/2,disc:="neutral"]                       # odd nb: middle bin = near-zero, colour grey
 H<-P[,.(med=median(revrate),n=.N),by=grp][order(grp)]
 lb<-rep("''",nb);lb[1]<-"atop('discordant','(strong)')";lb[nb]<-"atop('concordant','(strong)')";if(nb%%2==1)lb[(nb+1)/2]<-"rho[i]*rho[j]%~~%0"   # middle tertile labelled with the math expression, no 'neutral' word
 ggplot(P,aes(grpf,revrate))+geom_hline(yintercept=.5,linetype="dotted",color="grey65")+
  geom_violin(aes(fill=disc),color=NA,alpha=.22,scale="width",width=.85)+
  geom_jitter(aes(size=ns,color=disc),width=.17,alpha=.14,stroke=0)+
  geom_boxplot(width=.14,outlier.shape=NA,linewidth=.34,fill="white",color="grey30")+
  geom_point(data=H,aes(grp,med),shape=23,size=1.8,fill="white",color="black",stroke=.6)+
  scale_fill_manual(values=c(disc="#C0392B",conc="#2C6FBB",neutral="grey60"),guide="none")+scale_color_manual(values=c(disc="#C0392B",conc="#2C6FBB",neutral="grey55"),guide="none")+
  scale_size_continuous(range=c(.25,1.7),name="shared\ngenes")+
  scale_y_continuous(labels=function(x)x*100,limits=c(0,1))+scale_x_discrete(labels=parse(text=lb))+
  labs(x=expression("fertility coupling of disease pairs   ("*rho[i]%.%rho[j]*":  discordant "%->%" concordant)"),y="opposite-effect proportion (% shared genes)")+
  theme_classic(base_size=7)+theme(axis.text.x=element_text(size=5.7),axis.text.y=element_text(size=6),axis.title.y=element_text(size=6.5),axis.title.x=element_text(size=6,margin=margin(t=3)),legend.position="inside",legend.position.inside=c(.80,.96),legend.direction="horizontal",legend.title=element_text(size=5,vjust=.9),legend.key.size=unit(2.6,"mm"),legend.text=element_text(size=4.6),legend.margin=margin(0,0,0,0),plot.margin=margin(18,16,10,0))}
build<-function(sex,out,nb){R<-pipe(sprintf("ACAT_%s_gene_disease_all.csv",sex),sex);gA<-saddle_tri(R$P,R$rmap);gB<-binGlobal(copy(R$P),R$rmap,nb)
 m<-plot_grid(gA,gB,ncol=2,rel_widths=c(1.06,0.92),labels=c("A","B"),label_x=c(0.02,0.005),label_size=11,align="h",axis="t")
 ggsave(paste0(BD,out,".png"),m,width=182,height=88,units="mm",dpi=300,bg="white");ggsave(paste0(BD,out,".pdf"),m,width=182,height=88,units="mm",device=cairo_pdf,bg="white");cat(sex,nb,"->",out,"\n")}
build("female","Fig_intralocus_MAIN_AB_female_g3bin",3)
build("male","Fig_intralocus_SUPP_AB_male_g3bin",3)
