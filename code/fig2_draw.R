## Figure 2 draw — local sig-cell version. Left stacked (~25 dz) + CAD/SCZ/IBD zooms.
suppressMessages({library(data.table);library(tidyverse);library(ComplexHeatmap);library(circlize);library(grid)})
L<-readRDS("tmp/fig2_data.rds"); OUT<-"/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait"
disp<-L$disp; sc<-as.data.table(L$sc); dsum<-as.data.table(L$dsum)
## display-name cleanups
dsum[disease=="tsvColorectalcancer1",disp:="Colorectal_Cancer"][disp=="Lung_Carcinoma",disp:="Lung_Cancer"]
dsum[,disp:=gsub("_"," ",disp)]
## common tissue order (all tissues present), sorted
tiss<-sort(unique(sc$tissue))
col_fun<-colorRamp2(c(-10,0,15),c("#3B6FB6","white","#C84E4E"))
type_cols<-c("++"="#B2182B","--"="#2166AC","+-"="#C77CFF","-+"="#1B9E77")
class_cols<-c("ant"="#E64B35","syn"="#00A087")
MAXG<-35
## build one disease block: matrices genes x ALL tissues, rows split ant/syn
build_block<-function(dz,show_names=FALSE){
  gi<-disp[disease==dz]; scd<-sc[disease==dz & gene%in%gi$gene]
  if(!nrow(scd))return(NULL)
  # gene order: ant first, then by max|z|
  gr<-scd[,.(mz=max(abs(c(z_disease,z_fitness)),na.rm=TRUE)),by=gene]
  gr<-merge(gr,gi[,.(gene,conflict_type,conflict_class)],by="gene")
  gr<-gr[order(conflict_class!="ant",-mz)][1:min(.N,MAXG)]
  go<-gr$gene
  md<-dcast(scd[gene%in%go],gene~tissue,value.var="z_disease")
  mf<-dcast(scd[gene%in%go],gene~tissue,value.var="z_fitness")
  toM<-function(x){m<-as.matrix(x[,-1]);rownames(m)<-x$gene;out<-matrix(NA_real_,length(go),length(tiss),dimnames=list(go,tiss));cc<-intersect(colnames(m),tiss);out[rownames(m),cc]<-m[,cc];out[go,,drop=FALSE]}
  d<-dsum[disease==dz]
  list(md=toM(md),mf=toM(mf),ri=gr[match(go,gene)],
       label=sprintf("%s (%d/%d)",d$disp,d$n_ant,d$n_gene),show_names=show_names)
}
make_ht<-function(b,fs_title=7){
  rs<-factor(b$ri$conflict_class,levels=c("ant","syn"))
  ra<-rowAnnotation(Type=b$ri$conflict_type,Class=b$ri$conflict_class,
     col=list(Type=type_cols,Class=class_cols),show_annotation_name=FALSE,width=unit(4,"mm"),
     show_legend=FALSE)
  h1<-Heatmap(b$md,col=col_fun,na_col="grey95",cluster_rows=FALSE,cluster_columns=FALSE,
     row_split=rs,cluster_row_slices=FALSE,show_row_names=FALSE,show_column_names=FALSE,
     column_title=b$label,column_title_gp=gpar(fontsize=fs_title,fontface="bold"),
     row_title=NULL,rect_gp=gpar(col="white",lwd=0.2),border=TRUE,show_heatmap_legend=FALSE)
  h2<-Heatmap(b$mf,col=col_fun,na_col="grey95",cluster_rows=FALSE,cluster_columns=FALSE,
     row_split=rs,cluster_row_slices=FALSE,show_row_names=b$show_names,row_names_side="right",
     row_names_gp=gpar(fontsize=4.5),show_column_names=FALSE,column_title="male fertility",
     column_title_gp=gpar(fontsize=fs_title),row_title=NULL,rect_gp=gpar(col="white",lwd=0.2),
     border=TRUE,show_heatmap_legend=FALSE)
  ra+h1+h2
}
## ---- LEFT: select ~25 diseases (top ant_ratio per category, then top 25) ----
sel<-dsum[order(cat6,-ant_ratio,-n_ant)][,.SD[1:min(.N,5)],by=cat6][order(-ant_ratio)][1:min(.N,25)]
sel<-sel[order(factor(cat6,levels=c("Cardiometabolic","Immune","Neuropsychiatric","Cancer","Other")),-ant_ratio)]
cat("LEFT selected diseases:",nrow(sel),"\n"); print(sel[,.(disp,cat6,n_ant,n_gene)])
blocks<-lapply(sel$disease,build_block); blocks<-blocks[!sapply(blocks,is.null)]
htL<-make_ht(blocks[[1]]); for(i in 2:length(blocks)) htL<-htL %v% make_ht(blocks[[i]])
pdf(file.path(OUT,"Figure2_left_updated.pdf"),width=6,height=13)
draw(htL,padding=unit(c(3,3,3,3),"mm"),ht_gap=unit(1.2,"mm"))
draw(Legend(title="Z",col_fun=col_fun,direction="vertical"),x=unit(0.95,"npc"),y=unit(0.85,"npc"))
dev.off(); cat("saved Figure2_left_updated.pdf\n")
## ---- RIGHT: zoom CAD / SCZ / IBD with gene names ----
zoom<-c(cad="cad",scz="scz",ibd="ibd")
zb<-lapply(zoom,build_block,show_names=TRUE); zb<-zb[!sapply(zb,is.null)]
htZ<-make_ht(zb[[1]],fs_title=8); for(i in 2:length(zb)) htZ<-htZ %v% make_ht(zb[[i]],fs_title=8)
pdf(file.path(OUT,"Figure2_zoom_updated.pdf"),width=5,height=9)
draw(htZ,padding=unit(c(3,3,3,3),"mm"),ht_gap=unit(2,"mm"))
draw(Legend(title="Z",col_fun=col_fun,direction="vertical"),x=unit(0.93,"npc"),y=unit(0.9,"npc"))
draw(Legend(title="class",at=c("ant","syn"),labels=c("antagonistic","synergistic"),legend_gp=gpar(fill=class_cols)),x=unit(0.9,"npc"),y=unit(0.6,"npc"))
dev.off(); cat("saved Figure2_zoom_updated.pdf\n")
