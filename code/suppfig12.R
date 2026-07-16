## Supplementary Figures 1 & 2, IV-filtered consistent口径.
## Fig1: #significant disease associations per disease gene; all (orange) vs fertility-overlap (blue).
## Fig2: same, split by 4 disease categories, bins 1/2/3/>3.
suppressMessages({library(data.table);library(ggplot2)})
D<-"/Volumes/S840/04mypaper/lin/conflictresolution/dup_intralocus_analysis_2026-07-10"
OUT<-"/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait"
g<-fread(file.path(D,"SuppTable_IVfiltered_genes.csv"))
g<-g[disease_associated==TRUE]
g[,overlap:=(male_fertility|female_fertility)]          # blue subset = 267 dual
cat(sprintf("all disease genes=%d  fertility-overlap=%d\n",nrow(g),sum(g$overlap)))

## ---- category map (strip ' (ABBR)') ----
pd<-fread(file.path(D,"SuppTable_Fig1C_perdisease_IVfiltered.csv"))
strip<-function(x)trimws(gsub("\\s*\\([^)]*\\)","",x))
cat_of<-setNames(pd$Category,strip(pd$Disease))

## ---- per-gene long disease list ----
gl<-g[,.(dz=strip(unlist(strsplit(diseases,";")))),by=.(gene,overlap)]
gl[,cat:=cat_of[dz]]
cat("unmapped diseases:",gl[is.na(cat),uniqueN(dz)],"\n")

## ================= Supp Fig 1 =================
bins1<-function(n)factor(pmin(n,11),levels=1:11,labels=c(1:10,">10"))
cat(sprintf("dropping %d disease genes with n_diseases==0 (spurious, would form NA bar)\n",sum(g$n_diseases<1)))
d1<-g[n_diseases>=1,.(nd=n_diseases,grp=ifelse(overlap,"pleio","all"))]  # drop n_diseases==0 -> no NA bin
d1<-rbind(d1[,.(nd,grp)],d1[grp=="pleio"][,.(nd,grp="all")]) # 'all' = ALL disease genes (incl overlap), matching caption
d1[,bin:=bins1(nd)]
p1<-d1[,.N,by=.(grp,bin)][,prop:=N/sum(N),by=grp]
p1[,grp:=factor(grp,levels=c("all","pleio"),labels=c("all disease genes","+ fertility"))]
## stats for text
med_all<-median(g$n_diseases); med_pl<-median(g$n_diseases[g$overlap])
spec_all<-mean(g$n_diseases==1); spec_pl<-mean(g$n_diseases[g$overlap]==1); gt10_pl<-mean(g$n_diseases[g$overlap]>10)
cat(sprintf("MEDIAN: all=%d overlap=%d | disease-specific: all=%.1f%% overlap=%.1f%% | >10 overlap=%.1f%%\n",
  med_all,med_pl,100*spec_all,100*spec_pl,100*gt10_pl))

f1<-ggplot(p1,aes(bin,prop,fill=grp))+
  geom_col(position=position_dodge(preserve="single"),width=0.75)+
  scale_fill_manual(values=c("all disease genes"="#E8952C","+ fertility"="#3B6FB6"),name=NULL)+
  labs(x="Number of associated diseases",y="Proportion of genes")+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
    legend.position=c(0.78,0.85),legend.background=element_rect(fill=alpha("white",0.7),color=NA))
ggsave(file.path(OUT,"supp_fig1_pleiotropy.pdf"),f1,width=6,height=4,device=cairo_pdf)

## ================= Supp Fig 2 =================
cats<-c("Cancer","Cardiometabolic","Immune","Neuropsychiatric")
bins2<-function(n)factor(pmin(n,4),levels=1:4,labels=c("1","2","3",">3"))
## exclusive assignment: each gene -> its plurality category (sums to total disease genes, as in original)
gcc<-gl[,.N,by=.(gene,overlap,cat)]
prim<-gcc[order(gene,-N,cat)][,.SD[1],by=.(gene,overlap)]  # primary category + within-category disease count
gc<-prim[cat%in%cats,.(gene,overlap,cat,ncat=N)]           # per gene, #diseases in its (primary) category
lab<-c()
p2<-rbindlist(lapply(cats,function(cc){
  x<-gc[cat==cc]; n0<-uniqueN(x$gene); n1<-uniqueN(x[overlap==TRUE]$gene)
  lab[[cc]]<<-sprintf("%s\n(n=%d, %d)",cc,n0,n1)
  a<-x[,.(nd=ncat,grp="all")]; b<-x[overlap==TRUE,.(nd=ncat,grp="pleio")]
  d<-rbind(a,b); d[,bin:=bins2(nd)]
  r<-d[,.N,by=.(grp,bin)][,prop:=N/sum(N),by=grp][,cat:=cc]; r}))
p2[,facet:=factor(unlist(lab)[cat],levels=unlist(lab)[cats])]
p2[,grp:=factor(grp,levels=c("all","pleio"),labels=c("all disease genes","+ fertility"))]
cat("Supp Fig2 category sizes:\n");print(gc[,.(n0=uniqueN(gene),n1=uniqueN(gene[overlap])),by=cat])

f2<-ggplot(p2,aes(bin,prop,fill=grp))+
  geom_col(position=position_dodge(preserve="single"),width=0.75)+
  facet_wrap(~facet,nrow=2)+
  scale_fill_manual(values=c("all disease genes"="#E8952C","+ fertility"="#3B6FB6"),name=NULL)+
  labs(x="Number of associated diseases",y="Proportion of genes")+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
    strip.background=element_rect(fill="grey92"),legend.position="top")
ggsave(file.path(OUT,"supp_fig2_pleiotropy_bycat.pdf"),f2,width=7,height=5.2,device=cairo_pdf)
cat("saved supp_fig1_pleiotropy.pdf & supp_fig2_pleiotropy_bycat.pdf\n")
