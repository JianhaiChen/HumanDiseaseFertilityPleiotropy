## FusioS, aligned to the FusioMR standard: preclump/ (eQTL p<0.001) + >=5 IV per gene.
## One job = one endpoint x ALL tissues, so the FinnGen GWAS is parsed once, not 50 times.
## Usage: Rscript fusios_iv5.R <trait_gz> <pheno> [ntissue_chunk]
suppressMessages({library(data.table); library(parallel)})
CODE <- "${FUSIOMR_CODE}/fusiomr_sc_bk/code"
source("${BASEDIR}/fusios_finngen/code/functions_perchr.R")
source(file.path(CODE, "init_setup_seso.R"))
## SHARED cache: compile the Gibbs sampler once for the whole campaign,
## instead of once per array task (the old script keyed cacheDir on TASK_ID).
CC <- "${BASEDIR}/fusios_finngen/ccache_shared"
dir.create(CC, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(file.path(CODE, "gibbs_seso_uhp_only.cpp"), cacheDir = CC)

## Reproducibility: run_mr_seso is a Gibbs sampler with no internal seed, so two runs
## of the identical input differ by sampling noise. L'Ecuyer-CMRG gives mclapply
## reproducible per-worker streams (plain set.seed does NOT work across forks).
## NOTE: the stream depends on mc.cores, so keep --cpus-per-task fixed across the campaign.
RNGkind("L'Ecuyer-CMRG")
SEED <- 20260807L
set.seed(SEED)

a <- commandArgs(trailingOnly = TRUE)
trait_file <- a[1]; trait <- a[2]; ONE_TISSUE <- if (length(a) >= 3 && nzchar(a[3])) a[3] else NA_character_
NCORE <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
BASE <- "${BASEDIR}/fusios_finngen"
OUTB <- Sys.getenv("FUSIOS_OUTDIR", file.path(BASE, "tissue_sumtbl_iv5"))
MIN_IV <- 5                       # FusioMR: dt[, if (.N >= 5) .SD, by = gene]

TISSUES <- list.dirs(file.path(BASE, "preclump"), full.names = FALSE, recursive = FALSE)
TISSUES <- sort(TISSUES[nzchar(TISSUES)])
if (!is.na(ONE_TISSUE)) TISSUES <- intersect(TISSUES, ONE_TISSUE)

## ---- sex-limited endpoints must not borrow eQTLs from the other sex's tissue ----
FEMALE_ONLY_TISSUE <- c("Ovary","Uterus","Vagina","Fallopian_Tube",
                        "Cervix_Ectocervix","Cervix_Endocervix")
MALE_ONLY_TISSUE   <- c("Testis","Prostate")
sexmap <- tryCatch({
  sx <- fread("${BASEDIR}/endpoint_sex.tsv")
  setNames(sx$sex, sx$phenocode)
}, error = function(e) character(0))
esex <- if (length(sexmap) && trait %in% names(sexmap)) sexmap[[trait]] else "both"
if (esex == "female") {
  drop <- intersect(TISSUES, MALE_ONLY_TISSUE)
  if (length(drop)) cat(sprintf("  female-limited endpoint: dropping %s\n", paste(drop, collapse=",")))
  TISSUES <- setdiff(TISSUES, MALE_ONLY_TISSUE)
} else if (esex == "male") {
  drop <- intersect(TISSUES, FEMALE_ONLY_TISSUE)
  if (length(drop)) cat(sprintf("  male-limited endpoint: dropping %s\n", paste(drop, collapse=",")))
  TISSUES <- setdiff(TISSUES, FEMALE_ONLY_TISSUE)
}
if (!length(TISSUES)) { cat("SKIP no compatible tissue\n"); quit(status = 0) }
todo <- TISSUES[!vapply(TISSUES, function(t)
  file.exists(sprintf("%s/%s/fusios_%s_%s.txt", OUTB, t, trait, t)), logical(1))]
if (!length(todo)) { cat("SKIP all tissues done\n"); quit(status = 0) }
cat(sprintf("%s: %d/%d tissues to do, %d cores\n", trait, length(todo), length(TISSUES), NCORE))

## the expensive parse, done ONCE for all tissues
t0 <- Sys.time()
gwall <- fread(cmd = paste("gzip -dc", shQuote(trait_file)),
               select = c("#chrom","rsids","ref","alt","beta","sebeta"))
setnames(gwall, c("#chrom","rsids","ref","alt","beta","sebeta"),
         c("chrom","rsid","gw_oa","gw_ea","b_out","se_out"))
gwall[, chrom := suppressWarnings(as.integer(chrom))][, rsid := sub(",.*", "", rsid)]
gwall <- gwall[is.finite(b_out) & se_out > 0 & rsid != "" & !is.na(chrom)]
setkey(gwall, chrom)
cat(sprintf("  GWAS parsed %s rows in %.0fs\n", format(nrow(gwall), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

gene1 <- function(tmp) {
  if (nrow(tmp) < MIN_IV) return(NULL)          # skip BEFORE the Gibbs sampler
  bg <- max(1e-3, mean(tmp$beta^2) - mean(tmp$se^2))
  bt <- max(1e-5, mean(tmp$b_out^2) - mean(tmp$se_out^2))
  K <- nrow(tmp); ag <- at <- max(2, K/4)
  r <- tryCatch(run_mr_seso(tmp$gene_name[1], tmp$beta, tmp$se, tmp$b_out, tmp$se_out,
                            ag, bg*(ag-1), at, bt*(at-1)), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  c(r, b_gamma_prior_mean = bg, b_theta_prior_mean = bt)
}

for (tissue in todo) {
  ts <- Sys.time()
  set.seed(SEED)          # per-tissue reset: result must not depend on run order
  PC <- file.path(BASE, "preclump", tissue)
  outdir <- file.path(OUTB, tissue); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out <- sprintf("%s/fusios_%s_%s.txt", outdir, trait, tissue)
  if (file.exists(out)) next
  ACC <- list()
  for (chr in 1:22) {
    pf <- sprintf("%s/chr%d.txt", PC, chr); if (!file.exists(pf)) next
    pc <- tryCatch(fread(pf), error = function(e) NULL)
    if (is.null(pc) || !nrow(pc) || !"gene_name" %in% names(pc)) next
    ## drop under-instrumented genes up front -- ~6x fewer Gibbs calls
    pc <- pc[, if (.N >= MIN_IV) .SD, by = gene_name]
    if (!nrow(pc)) next
    gw <- gwall[.(chr), .(rsid, gw_ea, gw_oa, b_out, se_out), nomatch = 0L]
    d <- merge(pc, gw, by = "rsid")
    d[, align := fifelse(eqtl_ea == gw_ea & eqtl_oa == gw_oa, 1,
                  fifelse(eqtl_ea == gw_oa & eqtl_oa == gw_ea, -1, NA_real_))]
    d <- d[!is.na(align)][, b_out := b_out*align]
    d <- d[is.finite(beta) & is.finite(se) & is.finite(b_out) & is.finite(se_out)]
    if (!nrow(d)) next
    d <- d[, if (.N >= MIN_IV) .SD, by = gene_name]   # re-check after allele alignment
    if (!nrow(d)) next
    parts <- split(d, by = "gene_name", keep.by = TRUE)
    res <- if (NCORE > 1) mclapply(parts, gene1, mc.cores = NCORE) else lapply(parts, gene1)
    ## a crashed mclapply worker returns a try-error, NOT NULL.
    ## NOTE: run_mr_seso returns the gene NAME as its first element, so the result is a
    ## CHARACTER vector -- do not test is.numeric() here (that rejected every gene).
    bad <- vapply(res, function(x) is.null(x) || inherits(x, "try-error") || length(x) == 0,
                  logical(1))
    if (any(bad)) cat(sprintf("    chr%d: %d/%d genes failed\n", chr, sum(bad), length(res)))
    res <- do.call(rbind, res[!bad])
    if (is.null(res) || !length(res)) next
    res <- as.data.table(res); res[, chr := chr]; ACC[[length(ACC)+1]] <- res
  }
  final <- if (length(ACC)) rbindlist(ACC, fill = TRUE) else data.table()
  fwrite(final, out, sep = "\t")
  cat(sprintf("  %-42s rows %6d  %.0fs\n", tissue, nrow(final),
              as.numeric(difftime(Sys.time(), ts, units = "secs"))))
}
cat("DONE", trait, "\n")
