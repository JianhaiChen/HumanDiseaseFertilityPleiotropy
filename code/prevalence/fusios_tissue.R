## Sex-specific FusioS driver, generalized to ANY GTEx v10 tissue.
## Generalized from fusios_smoke.R (Muscle-hardcoded) -> tissue is now arg 1.
## Usage: Rscript fusios_tissue.R <tissue> <chr> <trait_summaryfile> <trait_short> <ngenes(0=all)>
suppressMessages({library(data.table); library(tidyverse); library(pbapply)})
CODE <- "/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code"
source(file.path(CODE, "functions.R"))
source(file.path(CODE, "init_setup_seso.R"))
Rcpp::sourceCpp(file.path(CODE, "gibbs_seso_uhp_only.cpp"),
  cacheDir = paste0("/scratch/jianhaichen/fusios_muscle/ccache/", Sys.getenv("SLURM_ARRAY_TASK_ID", "0")))

a <- commandArgs(trailingOnly = TRUE)
tissue <- a[1]; chr <- as.integer(a[2]); trait_file <- a[3]; trait <- a[4]; ng <- as.integer(a[5])
pcut <- 0.005; clump_kb <- 50; clump_r2 <- 0.01

eqf <- sprintf("/gpfs/data/linchen-lab/Bowei/dgtex/gtex_v10_rsid/eqtl_%s.v10.allpairs.chr%d.txt.gz", tissue, chr)
if (!file.exists(eqf)) { cat("NO EQTL", eqf, "\n"); quit(status=0) }
eq <- fread(eqf, select = c("gene_id","rsid","ref","alt","slope","slope_se","pval_nominal"))
eq <- eq[pval_nominal < 0.05 & is.finite(slope) & slope_se > 0]
eq[, gene_name := sub("\\..*","",gene_id)]
setnames(eq, c("ref","alt","slope","slope_se","pval_nominal"), c("eqtl_oa","eqtl_ea","beta","se","p"))

gw <- tryCatch(fread(trait_file), error=function(e) fread(cmd=paste("gzip -dc", shQuote(trait_file))))
gw <- gw[, .(rsid = variant_id, gw_ea = effect_allele, gw_oa = non_effect_allele,
             b_out = as.numeric(effect_size), se_out = as.numeric(standard_error))]
gw <- gw[is.finite(b_out) & is.finite(se_out) & se_out > 0]

d <- merge(eq[, .(gene_name, rsid, eqtl_ea, eqtl_oa, beta, se, p)], gw, by = "rsid")
d[, align := fifelse(eqtl_ea == gw_ea & eqtl_oa == gw_oa, 1,
              fifelse(eqtl_ea == gw_oa & eqtl_oa == gw_ea, -1, NA_real_))]
d <- d[!is.na(align)]; d[, b_out := b_out * align]
d <- d[is.finite(beta) & is.finite(se) & is.finite(b_out) & is.finite(se_out)]
cat("tissue", tissue, "chr", chr, "harmonized rows:", nrow(d), " genes:", uniqueN(d$gene_name), "\n")

genes <- unique(d$gene_name); if (ng > 0) genes <- head(genes, ng)
res <- do.call(rbind, pblapply(seq_along(genes), function(t)
  tryCatch(process_gene_bk_only(d, genes, t, pcut, clump_kb, clump_r2, NA, NA),
           error = function(e) NULL)))
res <- as.data.table(res)
outdir <- sprintf("/scratch/jianhaichen/fusios_muscle/sexspec_tissue/%s", tissue)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- sprintf("%s/fusios_%s_%s_chr%d.txt", outdir, trait, tissue, chr)
fwrite(res, out, sep = "\t")
cat("written:", out, "  rows:", nrow(res), "\n")
