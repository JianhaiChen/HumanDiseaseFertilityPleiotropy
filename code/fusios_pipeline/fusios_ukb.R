## FusioS UKB sexspec RECOMPUTE: article standard (P<=0.001 / 100kb / r2=0.1 / >=5 IV / top-50),
## per-chr bfile clump. Usage: Rscript fusios_ukb.R <trait_summaryfile> <trait> [chr_only]
suppressMessages({library(data.table); library(tidyverse); library(pbapply)})
CODE <- "/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code"
source("/scratch/jianhaichen/fusios_finngen/code/functions_perchr.R")   # bfile -> PERCHR_BFILE
source(file.path(CODE, "init_setup_seso.R"))
Rcpp::sourceCpp(file.path(CODE, "gibbs_seso_uhp_only.cpp"),
  cacheDir = paste0("/scratch/jianhaichen/fusios_muscle/ccache_v2/", Sys.getenv("SLURM_ARRAY_TASK_ID","0")))

a <- commandArgs(trailingOnly = TRUE)
trait_file <- a[1]; trait <- a[2]
chr_only <- if (length(a) >= 3 && a[3] != "") as.integer(a[3]) else NA_integer_
pcut <- 0.001; clump_kb <- 100; clump_r2 <- 0.1
outdir <- "/scratch/jianhaichen/fusios_muscle/sexspec_z_v2"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

## UKB summaryfile: variant_id effect_allele non_effect_allele effect_size standard_error
gwall <- tryCatch(fread(trait_file), error=function(e) fread(cmd=paste("gzip -dc", shQuote(trait_file))))
gwall <- gwall[, .(rsid = variant_id, gw_ea = effect_allele, gw_oa = non_effect_allele,
                   b_out = as.numeric(effect_size), se_out = as.numeric(standard_error))]
gwall <- gwall[is.finite(b_out) & is.finite(se_out) & se_out > 0 & rsid != ""]

chrset <- if (is.na(chr_only)) 1:22 else chr_only
for (chr in chrset) {
  PERCHR_BFILE <- sprintf("/scratch/jianhaichen/eur_perchr/chr%d", chr)
  out <- sprintf("%s/fusios_%s_Muscle_chr%d.txt", outdir, trait, chr)
  if (file.exists(out)) { cat("SKIP chr", chr, "\n"); next }
  eqf <- sprintf("/gpfs/data/linchen-lab/Bowei/dgtex/gtex_v10_rsid/eqtl_Muscle_Skeletal.v10.allpairs.chr%d.txt.gz", chr)
  eq <- fread(eqf, select = c("gene_id","rsid","ref","alt","slope","slope_se","pval_nominal"))
  eq <- eq[pval_nominal < 0.001 & is.finite(slope) & slope_se > 0]          # article: cis-eQTL P<=0.001
  eq[, gene_name := sub("\\..*","",gene_id)]
  setnames(eq, c("ref","alt","slope","slope_se","pval_nominal"), c("eqtl_oa","eqtl_ea","beta","se","p"))
  d <- merge(eq[, .(gene_name,rsid,eqtl_ea,eqtl_oa,beta,se,p)], gwall, by = "rsid")
  d[, align := fifelse(eqtl_ea==gw_ea & eqtl_oa==gw_oa, 1,
                fifelse(eqtl_ea==gw_oa & eqtl_oa==gw_ea, -1, NA_real_))]
  d <- d[!is.na(align)]; d[, b_out := b_out*align]
  d <- d[is.finite(beta) & is.finite(se) & is.finite(b_out) & is.finite(se_out)]
  ivn <- d[, .N, by = gene_name]; genes <- ivn[N >= 5, gene_name]           # article: >=5 IV
  if (length(genes) == 0) { fwrite(data.table(), out, sep="\t"); cat("EMPTY chr", chr, "\n"); next }
  d <- d[gene_name %in% genes]
  setorder(d, gene_name, p); d <- d[, head(.SD, 50L), by = gene_name]       # UKB: top-50 by eQTL p
  genes <- unique(d$gene_name)
  res <- do.call(rbind, pblapply(seq_along(genes), function(t)
    tryCatch(process_gene_bk_only(d, genes, t, pcut, clump_kb, clump_r2, NA, NA),
             error = function(e) NULL)))
  fwrite(as.data.table(res), out, sep = "\t")
  cat("chr", chr, "genes", length(genes), "rows", nrow(res), "\n")
}
cat("DONE", trait, "\n")
