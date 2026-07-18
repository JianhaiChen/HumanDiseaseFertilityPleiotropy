suppressMessages({library(data.table)});options(width=200)
strip<-function(x)sub("[.].*","",x)
setwd("/Volumes/X10Pro/data/synergistic/mr_res")
DICT<-fread("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait/disease_dictionary.csv")
co<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[,g:=strip(gene)];co<-co[chr%in%as.character(1:22)]
chrm<-setNames(as.integer(co$chr),co$g);posm<-setNames(as.numeric(co$pos),co$g)
d<-fread("ACAT_female_gene_disease_all.csv",select=c("gene","disease","b1_at_maxZ","padj_acat_1"))[,g:=strip(gene)][is.finite(b1_at_maxZ)]
d<-d[!(!is.na(chrm[g])&chrm[g]==6&posm[g]>25e6&posm[g]<34e6)]
d<-d[disease%in%DICT$acat_label];d[,ds:=sign(b1_at_maxZ)]
rf<-setNames(as.numeric(DICT$rhoGE_female),DICT$acat_label);rm<-setNames(as.numeric(DICT$rhoGE_male),DICT$acat_label)

pairtab<-function(rr){
  dz<-DICT$acat_label[is.finite(rr[DICT$acat_label])];dz<-intersect(dz,unique(d$disease))
  sig<-d[padj_acat_1<0.05&disease%in%dz,.(g,disease,ds)]
  Mx<-as.matrix(dcast(sig,g~disease,value.var="ds")[,-1]);rownames(Mx)<-dcast(sig,g~disease,value.var="ds")$g
  sc<-colnames(Mx);out<-list()
  for(i in 1:(length(dz)-1))for(j in (i+1):length(dz)){a<-dz[i];b<-dz[j];if(!(a%in%sc)||!(b%in%sc))next
    v<-Mx[,c(a,b)];v<-v[complete.cases(v),,drop=F];if(nrow(v)<15)next
    out[[length(out)+1]]<-data.table(a,b,ns=nrow(v),nrev=sum(v[,1]!=v[,2]),X=rr[a]*rr[b])}
  rbindlist(out)}

orbin<-function(sex,P){
  # bin pairs by |X| into quartiles; within each: OR = opp-effect odds(discordant X<0) / odds(concordant X>0)
  P[,absX:=abs(X)];P[,q:=cut(absX,quantile(absX,0:4/4),include.lowest=T,labels=c("Q1(weakest)","Q2","Q3","Q4(strongest)"))]
  res<-list()
  for(bb in levels(P$q)){Q<-P[q==bb];disc<-Q[X<0];conc<-Q[X>0]
    a11<-sum(disc$nrev);a12<-sum(disc$ns)-a11;a21<-sum(conc$nrev);a22<-sum(conc$ns)-a21
    or<-(a11*a22)/(a12*a21);se<-sqrt(1/a11+1/a12+1/a21+1/a22);ci<-exp(log(or)+c(-1.96,1.96)*se)
    res[[bb]]<-data.table(sex,bin=bb,npair=nrow(Q),
      medAbsX=signif(median(Q$absX),2),
      disc_opp=round(100*a11/(a11+a12),1),conc_opp=round(100*a21/(a21+a22),1),
      OR=round(or,2),lo=round(ci[1],2),hi=round(ci[2],2))}
  rbindlist(res)}
PF<-pairtab(rf);PM<-pairtab(rm)
cat("=== Female: OR by coupling-strength quartile ===\n");print(orbin("F",PF))
cat("\n=== Male ===\n");print(orbin("M",PM))
