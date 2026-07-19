## Aggregate per-gene z (=b_fusio/se_fusio) and compute rho_GE = cor(disease_z, father/mother_z)
## |z|<9 cap, >=30 shared genes. (MHC not excluded yet — preliminary.)
suppressMessages(library(data.table))
FG<-"/scratch/jianhaichen/fusios_finngen/sumtbl"
FM<-"/scratch/jianhaichen/fusios_muscle/sexspec_z_v2"
getz<-function(dir,trait){
  fs<-sprintf("%s/fusios_%s_Muscle_chr%d.txt",dir,trait,1:22); fs<-fs[file.exists(fs)]
  x<-rbindlist(lapply(fs,function(f){tryCatch(fread(f),error=function(e)NULL)}),fill=TRUE)
  if(is.null(x)||nrow(x)==0||!"b_fusio"%in%names(x)) return(NULL)
  x<-x[is.finite(b_fusio)&is.finite(se_fusio)&se_fusio>0]; x[,z:=b_fusio/se_fusio]; x<-x[abs(z)<9]
  unique(x[,.(gene,z)],by="gene")
}
fa<-getz(FM,"father"); mo<-getz(FM,"mother")
setnames(fa,"z","zf"); setnames(mo,"z","zm")
cat("father genes:",nrow(fa)," mother genes:",nrow(mo),"\n")
allf<-list.files(FG,pattern="_Muscle_chr[0-9]+.txt")
traits<-unique(sub("fusios_(.*)_Muscle_chr[0-9]+\\.txt","\\1",allf))
done<-traits[sapply(traits,function(t) sum(file.exists(sprintf("%s/fusios_%s_Muscle_chr%d.txt",FG,t,1:22)))==22)]
cat("completed endpoints:",length(done),"\n")
res<-rbindlist(lapply(done,function(t){
  z<-getz(FG,t); if(is.null(z)||nrow(z)<30) return(NULL)
  df<-merge(z,fa,by="gene"); dm<-merge(z,mo,by="gene")
  data.table(phenocode=t,
    rho_father=if(nrow(df)>=30) cor(df$z,df$zf) else NA_real_, ng_f=nrow(df),
    rho_mother=if(nrow(dm)>=30) cor(dm$z,dm$zm) else NA_real_, ng_m=nrow(dm))
}),fill=TRUE)
fwrite(res,"/scratch/jianhaichen/fusios_finngen/rhoGE_partial.tsv",sep="\t")
cat("written",nrow(res),"endpoints with rho_GE\n")
