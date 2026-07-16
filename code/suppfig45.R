## Supp Fig 4 (overall) & Supp Fig 5 (by disease category): 4 effect-direction classes,
## proportions within each fertility trait, counts above bars. IV-filtered kept dual set.
suppressMessages({library(data.table);library(ggplot2)})
MR<-"/Volumes/X10Pro/data/synergistic/mr_res"; OUT<-"/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait"
k<-readRDS("/Users/cjh/.claude/jobs/bb3be850/tmp/antag_k.rds")
info<-fread(file.path(MR,"disease_gwas_source_info.csv"),select=c("disease","Category"))
mapcat<-function(c)fifelse(c=="neuro_psych","Neuropsychiatric",fifelse(c=="immune","Immune",
  fifelse(c%in%c("cardiovasc","metabolic"),"Cardiometabolic",fifelse(c=="cancer","Cancer",
  fifelse(c=="reproductive","Reproductive","Other")))))
info[,Category:=mapcat(Category)]; k<-merge(k,info,by="disease",all.x=TRUE)
k[,Sex:=fifelse(fitness=="male","Male","Female")]
code<-c("d+f+"="++","d-f-"="--","d-f+"="-+","d+f-"="+-")
k[,cc:=code[klass]]
lev<-c("++","--","+-","-+")
lab4<-c("++"="d+f+ (antag.)","--"="d-f- (antag.)","+-"="d+f- (syn.)","-+"="d-f+ (syn.)")
cols4<-c("++"="#C0392B","--"="#E67E22","+-"="#2C5985","-+"="#5DADE2")
k[,cc:=factor(cc,levels=lev)]
th<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
  legend.position="top",strip.background=element_rect(fill="grey92"))

## ---- Supp Fig 4: overall, within each sex ----
d4<-k[,.N,by=.(Sex,cc)][,prop:=N/sum(N),by=Sex]
f4<-ggplot(d4,aes(cc,prop,fill=cc))+geom_col(width=0.75)+
  geom_text(aes(label=N),vjust=-0.3,size=3)+facet_wrap(~Sex)+
  scale_fill_manual(values=cols4,labels=lab4,name=NULL)+
  scale_y_continuous(expand=expansion(mult=c(0,0.12)),labels=scales::percent)+
  labs(x="Effect-direction class",y="Proportion of pleiotropic gene-disease pairs")+th
ggsave(file.path(OUT,"supp_fig3_4class_overall.pdf"),f4,width=6,height=4,device=cairo_pdf)

## ---- Supp Fig 5: within each sex x disease category ----
cats<-c("Cancer","Cardiometabolic","Immune","Neuropsychiatric","Other","Reproductive")
d5<-k[Category%in%cats,.N,by=.(Sex,Category,cc)][,prop:=N/sum(N),by=.(Sex,Category)]
d5[,Category:=factor(Category,levels=cats)]
f5<-ggplot(d5,aes(cc,prop,fill=cc))+geom_col(width=0.78)+
  geom_text(aes(label=N),vjust=-0.3,size=2.4)+
  facet_grid(Sex~Category)+
  scale_fill_manual(values=cols4,labels=lab4,name=NULL)+
  scale_y_continuous(expand=expansion(mult=c(0,0.15)),labels=scales::percent)+
  labs(x="Effect-direction class",y="Proportion of pleiotropic gene-disease pairs")+
  th+theme(axis.text.x=element_text(size=7,angle=45,hjust=1))
ggsave(file.path(OUT,"supp_fig4_4class_bycat.pdf"),f5,width=10,height=5,device=cairo_pdf)
cat("saved supp_fig3_4class_overall.pdf & supp_fig4_4class_bycat.pdf\n")
cat("\nSupp Fig4 within-sex proportions:\n");print(dcast(d4,Sex~cc,value.var="prop")[,lapply(.SD,function(x)if(is.numeric(x))round(100*x,1) else x)])
