## cis-eQTLs p<0.001, LD-clumped in R (100 kb, r2<0.1). No cap per gene.
## Usage: Rscript preclump_R_full.R <tissue> [chr]
## Pre-clump via in-R LD (BEDMatrix). Usage: Rscript preclump_R.R <tissue> [chr]
suppressMessages({library(data.table); library(BEDMatrix)})
a <- commandArgs(trailingOnly=TRUE); tissue <- a[1]; chr_only <- if(length(a)>=2 && a[2]!="") as.integer(a[2]) else NA_integer_
KB <- 100000; R2T <- 0.1
ANNOT <- "${BASEDIR}/rerun_conflict/gtex_v10_annot"
REF <- "${BASEDIR}/finngen_r13/finngen_R13_AB1_ACTINOMYCOSIS.gz"
outdir <- sprintf("${WORKDIR}/preclump_full/%s", tissue); dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
LBF <- sprintf("/dev/shm/pcRbf_%s_%s", Sys.getenv("SLURM_JOB_ID","x"), Sys.getenv("SLURM_ARRAY_TASK_ID","0")); dir.create(LBF, showWarnings=FALSE, recursive=TRUE); on.exit(unlink(LBF, recursive=TRUE), add=TRUE)
refgw <- fread(cmd=paste("gzip -dc", shQuote(REF)), select=c("#chrom","rsids","ref","alt"))
setnames(refgw, c("#chrom","rsids","ref","alt"), c("chrom","rsid","gw_oa","gw_ea"))
refgw[, chrom := suppressWarnings(as.integer(chrom))][, rsid := sub(",.*","",rsid)]; refgw <- refgw[!is.na(chrom) & rsid!=""]
clumpR <- function(rs, pp, G, posv){
  ok <- rs[rs %in% colnames(G)]; if(length(ok)<1) return(character(0))
  ok <- ok[order(pp[match(ok,rs)])]; keep <- character(0); rem <- ok
  while(length(rem)>0){ top <- rem[1]; keep <- c(keep,top); rem <- rem[-1]; if(length(rem)==0) break
    win <- abs(posv[rem]-posv[top]) <= KB
    if(any(win)){ r2 <- as.numeric(suppressWarnings(cor(G[,top], G[,rem[win],drop=FALSE], use="pairwise.complete.obs")))^2
      r2[is.na(r2)] <- 0; drop <- rem[win][r2>R2T]; rem <- rem[!rem %in% drop] } }
  keep }
chrset <- if(is.na(chr_only)) 1:22 else chr_only
for(chr in chrset){
  out <- sprintf("%s/chr%d.txt", outdir, chr); if(file.exists(out) && is.na(chr_only)) next
  bpref <- sprintf("%s/chr%d", LBF, chr); if(!file.exists(paste0(bpref,".bed"))) invisible(file.copy(Sys.glob(sprintf("${BASEDIR}/eur_perchr/chr%d.*",chr)), LBF))
  eqf <- sprintf("%s/GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-%s.v10.allpairs.chr%d.annot.txt.gz", ANNOT, tissue, chr)
  if(!file.exists(eqf)){ cat("NO EQTL",tissue,chr,"\n"); next }
  eq <- fread(eqf, select=c("gene_id_noversion","rsid","effect_allele","ref_allele","slope","slope_se","pval"))[pval<0.001 & is.finite(slope) & slope_se>0]
  setnames(eq, c("gene_id_noversion","effect_allele","ref_allele","slope","slope_se","pval"), c("gene_name","eqtl_ea","eqtl_oa","beta","se","p"))
  gw <- refgw[chrom==chr]
  d <- merge(eq[,.(gene_name,rsid,eqtl_ea,eqtl_oa,beta,se,p)], gw[,.(rsid,gw_ea,gw_oa)], by="rsid")
  d[, align := fifelse(eqtl_ea==gw_ea & eqtl_oa==gw_oa, 1, fifelse(eqtl_ea==gw_oa & eqtl_oa==gw_ea, -1, NA_real_))]; d <- d[!is.na(align)]
  ivn <- d[,.N,by=gene_name]; genes <- ivn[N>=5, gene_name]; if(length(genes)==0){ fwrite(data.table(), out, sep="\t"); next }
  d <- d[gene_name %in% genes]; setorder(d, gene_name, p); genes <- unique(d$gene_name)
  bim <- fread(paste0(bpref,".bim")); setnames(bim, 1:6, c("c","rsid","cm","pos","a1","a2"))
  bm <- BEDMatrix(bpref)
  ivu <- unique(d$rsid); mi <- match(ivu, bim$rsid); ok <- !is.na(mi); ivu <- ivu[ok]; mi <- mi[ok]
  G <- bm[, mi, drop=FALSE]; colnames(G) <- ivu; posv <- setNames(bim$pos[mi], ivu)
  res <- rbindlist(lapply(genes, function(g){ tmp <- d[gene_name==g]; kp <- clumpR(tmp$rsid, tmp$p, G, posv)
    if(length(kp)<1) return(NULL); tmp[rsid %in% kp, .(gene_name,rsid,eqtl_ea,eqtl_oa,beta,se,p)] }), fill=TRUE)
  fwrite(res, out, sep="\t"); cat("chr",chr,"genes",length(genes),"clumped_ivs",nrow(res),"\n")
}
cat("PRECLUMP_R DONE",tissue,"\n")
