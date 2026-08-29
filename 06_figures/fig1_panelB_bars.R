# Figure 1c — IV-filtered consistent: disease genes per-disease IV-filtered; fertility = 182/169 (IV-filtered consistent).
# names+colors matched to fig1c_rerun_mirror_97.pdf.
suppressMessages({library(data.table);library(ggplot2);library(scales)}); setwd("${MR}")
SC<-"${PROJ}/figures_final"
RECOMB<-"${MR}/recomb"
cm_gap<-0.5;sh_n<-5L;sh_frac<-0.5;cor_same<-0.8;cor_flip<--0.8
co<-fread("gene_coords_v19.tsv",header=FALSE,col.names=c("gene","chr","pos"))[chr%in%as.character(1:22)][,.(gene,chr=as.integer(chr),pos=as.numeric(pos))]
cmf<-vector("list",22);for(ch in 1:22){f<-sprintf("%s/plink.chr%d.GRCh37.map",RECOMB,ch);if(file.exists(f)){m<-fread(f)[,.(bp=V4,cm=V3)][order(bp)];cmf[[ch]]<-approxfun(m$bp,m$cm,rule=2)}}
co[,cm:=as.numeric(NA)];for(ch in 1:22)if(!is.null(cmf[[ch]]))co[chr==ch,cm:=cmf[[ch]](pos)];co[,in_mhc:=chr==6&pos>25e6&pos<34e6]
pr<-fread("session_2026-07-09_iv_filter/pairs_fdr05.csv",select=c("geneA","geneB","nA","nB","shared_n","eqtl_cor"))
pr[,`:=`(a=pmin(geneA,geneB),b=pmax(geneA,geneB))];pl<-pr[order(-shared_n)][,.SD[1],by=.(a,b)][,.(gene_a=a,gene_b=b,shared_n,nA,nB,eqtl_cor)]
pl[,shared_frac:=shared_n/pmin(nA,nB)];pl[,high:=shared_n>=sh_n&!is.na(shared_frac)&shared_frac>=sh_frac];pl[,flip:=high&!is.na(eqtl_cor)&eqtl_cor<=cor_flip];pl[,redun:=high&!is.na(eqtl_cor)&eqtl_cor>=cor_same];setkey(pl,gene_a,gene_b)
fc<-function(nodes,edges){nodes<-sort(unique(nodes));par<-setNames(nodes,nodes);fr<-function(x){while(par[x]!=x)x<-par[x];x};if(nrow(edges))for(i in seq_len(nrow(edges))){ra<-fr(edges$gene_a[i]);rb<-fr(edges$gene_b[i]);if(ra!=rb)par[ra]<-rb};data.table(gene=nodes,component=sapply(nodes,fr))}
filt<-function(dt){dt<-merge(dt,co,by="gene",all.x=TRUE);dt[,miss:=is.na(chr)|is.na(pos)|is.na(cm)];out<-list()
 for(st in unique(dt$stratum)){x0<-dt[stratum==st];out[[length(out)+1]]<-x0[miss==TRUE][,act:="keep_no_coord"];out[[length(out)+1]]<-x0[miss==FALSE&in_mhc==TRUE][,act:="drop_mhc"];x<-x0[miss==FALSE&in_mhc==FALSE];if(!nrow(x))next
  setorder(x,chr,cm,pos,gene);x[,blk:=cumsum(is.na(shift(chr))|chr!=shift(chr)|(cm-shift(cm))>cm_gap)]
  for(b in unique(x$blk)){xb<-x[blk==b];g<-sort(unique(xb$gene));xb[,act:="keep"];signs<-sort(unique(na.omit(xb$sign)));patt<-if(length(signs)==1)"same" else if(length(signs)>1)"mixed" else "miss"
   if(length(g)==1){xb[,act:="keep_singleton"]}else{gp<-CJ(gene_a=g,gene_b=g)[gene_a<gene_b];bp<-pl[gp,on=.(gene_a,gene_b)];bp[is.na(high),`:=`(high=FALSE,flip=FALSE,redun=FALSE)];sg<-setNames(xb$sign,xb$gene);bp[,`:=`(sa=sg[gene_a],sb=sg[gene_b])];bp[,`:=`(opp=!is.na(sa)&!is.na(sb)&sa!=sb,same=!is.na(sa)&!is.na(sb)&sa==sb)]
    if(patt!="mixed"){re<-bp[redun==TRUE&same==TRUE,.(gene_a,gene_b)];if(nrow(re)){cp<-fc(g,re);xb<-merge(xb,cp,by="gene",all.x=TRUE);for(cc in cp[,.N,by=component][N>1,component]){mem<-xb[component==cc];rp<-mem[order(-strength,gene)]$gene[1];xb[component==cc&gene==rp,act:="keep_rep"];xb[component==cc&gene!=rp,act:="drop_redundant_same_sign"]};xb[act=="keep",act:="keep_same"];xb[,component:=NULL]}else xb[,act:="keep_same"]}
    else{conf<-bp[flip==TRUE&opp==TRUE,.(gene_a,gene_b)];poss<-bp[high==TRUE&opp==TRUE&flip!=TRUE,.(gene_a,gene_b)];de<-unique(rbind(conf[,lab:="c"],poss[,lab:="p"]));if(nrow(de)){cp<-fc(g,de[,.(gene_a,gene_b)]);xb<-merge(xb,cp,by="gene",all.x=TRUE);ec<-merge(de,cp[,.(gene_a=gene,component)],by="gene_a",all.x=TRUE);cl<-ec[,.(hc=any(lab=="c"),hp=any(lab=="p")),by=component];xb[component%in%cl[hc==TRUE,component],act:="drop_confirmed_mixed_artifact"];xb[component%in%cl[hc==FALSE&hp==TRUE,component],act:="drop_possible_mixed_artifact"];xb[act=="keep",act:="keep_mixed"];xb[,component:=NULL]}else xb[,act:="keep_mixed"]}}
   out[[length(out)+1]]<-xb}}
 A<-rbindlist(out,use.names=TRUE,fill=TRUE);A[,keep:=!act%in%c("drop_mhc","drop_redundant_same_sign","drop_confirmed_mixed_artifact","drop_possible_mixed_artifact")];A}
## fertility set 182/169
iv<-fread("bonf050_broadmixed_recomb_block_shared_iv_filter_combos_with_annotation.csv");ivk<-iv[keep_after_filter==TRUE];allivg<-unique(iv$gene)
b2<-fread("gene_level_b2_fitness_acat_summary.dropbox.csv",select=c("gene","fitness","p_bonf_global","meta_beta_ivw","meta_z_ivw","top_beta","top_se"))[!is.na(fitness)&p_bonf_global<0.05]
fg<-list()
for(s in c("male","female")){D<-unique(ivk[fitness==s]$gene);fspec<-setdiff(unique(b2[fitness==s]$gene),allivg);fb<-b2[fitness==s&gene%in%fspec];fb[,`:=`(stratum="F",sign=fifelse(meta_beta_ivw>0,"+","-"),strength=fifelse(!is.na(meta_z_ivw),abs(meta_z_ivw),abs(top_beta/top_se)))];fk<-if(nrow(fb))unique(filt(fb[,.(gene,stratum,sign,strength)])[keep==TRUE]$gene) else character();fg[[s]]<-union(D,fk)}
## disease per-disease IV-filtered
b1<-fread("gene_level_b1_disease_acat_summary.dropbox.csv",select=c("gene","disease","p_bonf_global","meta_beta_ivw","meta_z_ivw","top_beta","top_se"))[p_bonf_global<0.05]
b1[,`:=`(stratum=disease,sign=fifelse(meta_beta_ivw>0,"+","-"),strength=fifelse(!is.na(meta_z_ivw),abs(meta_z_ivw),abs(top_beta/top_se)))]
dk<-filt(b1[,.(gene,stratum,sign,strength)])[keep==TRUE,.(gene,disease=stratum)]
## names+colors (match original)
info<-read.csv("disease_gwas_source_info.csv",stringsAsFactors=FALSE)[,c("disease","Disease","Category")]
cn<-function(s){s<-gsub("\\s*\\([^)]*\\)","",s);gsub("[ /]+","_",trimws(s))}
mc<-function(c)ifelse(c=="neuro_psych","Neuropsychiatric",ifelse(c=="immune","Immune",ifelse(c%in%c("cardiovasc","metabolic"),"Cardiometabolic",ifelse(c=="cancer","Cancer",ifelse(c=="reproductive","Reproductive","Other")))))
disp<-setNames(cn(info$Disease),info$disease);cat0<-setNames(mc(info$Category),info$disease)
origset<-c("Bladder_Cancer","Prostate_Cancer","Leukemia","Colorectal_Cancer","Head_Neck_Cancer","Skin_Cancer","Brain_Cancer","Renal_Cancer","Multiple_Myeloma","Liver_Cancer","Breast_Cancer","Thyroid_Cancer","Lung_Cancer","Endometrial_Cancer","Hodgkins_Lymphoma","Esophageal_Cancer","Non-Hodgkins_Lymphoma","Ovarian_Cancer","Melanoma","Gastric_Cancer","Pancreatic_Cancer","Aortic_Valve_Stenosis","Hypercholesterolemia","Obesity","Gout","Coronary_Artery_Disease","Hypertension","Atrial_Fibrillation","Myocardial_Infarction","Varicose_Veins","Abnormal_Thrombosis","Deep_Vein_Thrombosis","Blood_Clot_Lung","Type_2_Diabetes","Hypertrophic_Cardiomyopathy","Heart_Failure","Stroke","Raynauds_Disease","Gestational_Diabetes","VitaminD_Deficiency","Sepsis","Sarcoidosis","Graves_Disease","Acne","Rheumatoid_Arthritis","Asthma_Pneumonia","Hypothyroidism","Rosacea","Psoriasis","Hay_Fever","Type_1_Diabetes","Inflammatory_Bowel_Disease","Ulcerative_Colitis","Alopecia_Areata","Uveitis","Systemic_Lupus_Erythematosus","Tuberculosis","Crohns_Disease","Atopic_Dermatitis_Eczema","Myasthenia_Gravis","Primary_Biliary_Cirrhosis","Schizophrenia","Restless_Legs_Syndrome","Panic_Disorder","Bipolar_Disorder","Insomnia","Frontotemporal_Dementia","Alzheimers_Disease","Obstructive_Sleep_Apnea","Tourette_Syndrome","Borderline_Personality_Disorder","Lewy_Body_Dementia","Parkinsons_Disease","ADHD","PTSD","Autism_Spectrum_Disorder","Major_Depressive_Disorder","Iron_Deficiency_Anemia","Gallstones","Macular_Degeneration","COPD","OpenAngle_Glaucoma","Nephrotic_Syndrome","Osteoporosis","Hearing_Loss","Cataract","Osteoarthritis","Chronic_Kidney_Disease","Kidney_Stones","Idiopathic_Pulmonary_Fibrosis","High_Myopia","Benign_Prostatic_Hyperplasia","Polycystic_Ovary_Syndrome","Uterine_Fibroids","Endometriosis","Preeclampsia","Testicular_Disease")
nz<-function(s)tolower(gsub("[^a-z0-9]","",s));obk<-setNames(origset,nz(origset))
raw<-sort(unique(dk$disease))
dispmap<-setNames(vapply(raw,function(r){d<-unname(disp[r]);h<-unname(obk[nz(d)]);if(!is.na(h))h else if(!is.na(d))d else r},character(1)),raw)
dispmap["tsvColorectalcancer1"]<-"Colorectal_Cancer";dispmap["lungcarcinoma"]<-"Lung_Cancer";dispmap["Blood_Clot_Lung"]<-"Blood_Clot_Lung"
ab<-c("Hypertrophic_Cardiomyopathy"="Hypertrophic_Cardiomyopathy (HCM)","Inflammatory_Bowel_Disease"="Inflammatory_Bowel_Disease (IBD)","Systemic_Lupus_Erythematosus"="Systemic_Lupus_Erythematosus (SLE)","Restless_Legs_Syndrome"="Restless_Legs_Syndrome (RLS)","Frontotemporal_Dementia"="Frontotemporal Dementia (FTD)","Borderline_Personality_Disorder"="Borderline_Personality_Disorder (BPD)","Autism_Spectrum_Disorder"="Autism Spectrum Disorder (ASD)","Major_Depressive_Disorder"="Major_Depressive_Disorder (MDD)","Idiopathic_Pulmonary_Fibrosis"="Idiopathic_Pulmonary_Fibrosis (IPF)","Benign_Prostatic_Hyperplasia"="Benign_Prostatic_Hyperplasia (BPH)","Polycystic_Ovary_Syndrome"="Polycystic_Ovary_Syndrome (PCOS)")
dispmap<-setNames(ifelse(dispmap%in%names(ab),unname(ab[dispmap]),dispmap),names(dispmap))
catmap<-cat0[raw];catmap[is.na(catmap)]<-"Cancer";names(catmap)<-raw
fem<-c("EFO_0000305Breastcancer.h.tsv.gz1","endometrialcancerwu","Endometriosis","Gestational_diabetes","ovariancancerwu","Polycystic_ovary_syndrome","Preeclampsia","Uterine_Fibroids")
mal<-c("GCST90274714Prostatecancer","Hyperplasiaofprostate","Testicular_disease")
cat_order<-c("Reproductive","Other","Neuropsychiatric","Immune","Cardiometabolic","Cancer")
one<-function(sex){dds<-setdiff(raw,if(sex=="male")fem else mal);do.call(rbind,lapply(dds,function(dd){dg<-dk[disease==dd]$gene;nd<-length(dg);ns<-length(intersect(dg,fg[[sex]]));data.frame(disease=dd,Disease=unname(dispmap[dd]),Category=unname(catmap[dd]),nd=nd,ns=ns,prop_d=if(nd>0)ns/nd else NA,sex=tools::toTitleCase(sex),stringsAsFactors=FALSE)}))}
res<-rbind(one("male"),one("female"))
cat(sprintf("male median=%.1f%% female median=%.1f%%\n",100*median(res$prop_d[res$sex=="Male"],na.rm=T),100*median(res$prop_d[res$sex=="Female"],na.rm=T)))
agg<-aggregate(prop_d~Disease+Category,res,max,na.rm=TRUE);agg$Category<-factor(agg$Category,levels=cat_order);agg<-agg[order(agg$Category,-agg$prop_d,agg$Disease),];ord<-agg$Disease
res$Disease<-factor(res$Disease,levels=ord)
dm<-res[res$sex=="Male"&!res$disease%in%fem,];dm$value<--dm$prop_d;dm$lbl<-sprintf("%d/%d",dm$ns,dm$nd)
df_<-res[res$sex=="Female"&!res$disease%in%mal,];df_$value<-df_$prop_d;df_$lbl<-sprintf("%d/%d",df_$ns,df_$nd)
d_all<-rbind(dm,df_);d_all$Category<-factor(d_all$Category,levels=cat_order);df_axis<-data.frame(Disease=factor(ord,levels=ord),y=0);pad<-0.012
saveRDS(list(dk=dk,fg=fg,res=res),file.path(SC,"fig1_panelB_data.rds"))
cat("raw:",length(raw)," dispmap NA:",sum(is.na(dispmap))," dup disp:",sum(duplicated(dispmap))," res rows:",nrow(res)," res Disease NA:",sum(is.na(as.character(res$Disease)))," ord len:",length(ord),"\n")
cat("dup display names:",paste(dispmap[duplicated(dispmap)|duplicated(dispmap,fromLast=T)],collapse=","),"\n")
pp<-ggplot(df_axis,aes(Disease,y))+geom_col(data=d_all,aes(y=value,fill=Category,alpha=sex),width=0.5,linewidth=0.4)+
  geom_text(data=d_all,aes(y=value+ifelse(sex=="Male",-pad,pad),label=lbl,hjust=ifelse(sex=="Male",1,0)),size=2)+geom_hline(yintercept=0,col="gray40")+
  scale_alpha_manual(values=c(Male=0.55,Female=0.9),guide="none")+
  scale_fill_manual(values=c(Reproductive="#E86DE8",Other="#A6A6A6",Neuropsychiatric="#7DA0F9",Immune="#3CCB8B",Cardiometabolic="#E89B2C",Cancer="#F47C7C"))+
  guides(fill="none")+
  scale_y_continuous(limits=c(-0.5,0.5),labels=function(x)percent(abs(x),accuracy=1),breaks=breaks_pretty(n=6))+scale_x_discrete(drop=FALSE,limits=ord)+coord_flip(clip="off")+
  labs(x=NULL,y="% of disease genes with fertility effects")+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title=element_text(size=8),axis.text=element_text(size=6),legend.title=element_text(size=8,face="bold"),legend.text=element_text(size=7),legend.position="right")+
  annotate("text",x=c(length(ord)-6,length(ord)-9.5),y=c(-0.47,0.47),label=c("Male","Female"),hjust=c(1,0),angle=90,size=3)
ggsave(file.path(SC,"Fig1_panelB_bars.pdf"),pp,width=3.35,height=8.6,device=cairo_pdf);ggsave(file.path(SC,"Fig1_panelB_bars.png"),pp,width=3.35,height=8.6,dpi=150);cat("saved\n")
