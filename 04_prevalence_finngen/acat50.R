## 50-tissue gene-level ACAT for FinnGen endpoints, with sex-tissue compatibility.
## ACAT formula follows the canonical gene_level_acat_summary.R (tan-sum / atan).
## Usage: Rscript acat50.R <chunk_index> <n_chunks>
suppressMessages(library(data.table))
BASE <- "${BASEDIR}/fusios_finngen"
TS   <- file.path(BASE, "tissue_sumtbl_iv5")
OUT  <- "${BASEDIR}/acat50"; dir.create(OUT, showWarnings = FALSE)

## Sex-limited tissue lists MUST match fusios_iv5.R exactly. acat26.R listed only
## Ovary/Uterus and Testis, so Prostate could feed gynaecological endpoints and
## Vagina could feed male ones. In practice fusios_iv5.R already writes no file
## for those combinations, but the mismatch corrupted n_tis_avail and would bite
## the moment anyone reran the tissue step without the sex rule.
FEMALE_ONLY <- c("Ovary","Uterus","Vagina","Fallopian_Tube",
                 "Cervix_Ectocervix","Cervix_Endocervix")
MALE_ONLY   <- c("Testis","Prostate")

## Expected endpoint count per tissue, derived not hardcoded: endpoint_sex.tsv
## lists more sex-limited endpoints than endpoints.tsv actually contains.
eps_all <- fread("${BASEDIR}/endpoints.tsv", header = FALSE)$V1
sx <- fread("${BASEDIR}/endpoint_sex.tsv")
sexmap <- setNames(sx$sex, sx$phenocode)
n_fem_ep <- sum(eps_all %in% sx$phenocode[sx$sex == "female"])
n_mal_ep <- sum(eps_all %in% sx$phenocode[sx$sex == "male"])
expected <- function(t) {
  if (t %in% FEMALE_ONLY) length(eps_all) - n_mal_ep
  else if (t %in% MALE_ONLY) length(eps_all) - n_fem_ep
  else length(eps_all)
}

## A tissue joins the ACAT only if its run actually finished. A half-finished
## tissue would silently shrink n_tissue for exactly the endpoints it is missing,
## which biases the Cauchy combination rather than just adding noise.
MIN_FRAC <- 0.98
TISSUES <- list.dirs(TS, full.names = FALSE, recursive = FALSE)
TISSUES <- TISSUES[nzchar(TISSUES)]
ok <- vapply(TISSUES, function(t)
  length(list.files(file.path(TS, t), pattern = "\\.txt$")) >= MIN_FRAC * expected(t),
  logical(1))
if (any(!ok)) cat("EXCLUDED (incomplete):", paste(TISSUES[!ok], collapse = ", "), "\n")
TISSUES <- TISSUES[ok]
cat(sprintf("tissues used: %d\n", length(TISSUES)))

## Drop a tissue's contribution for a gene when that tissue had too few
## instruments to be trustworthy. The pipeline already enforces >=5, so this only
## bites if MIN_IV is ever loosened upstream.
MIN_IV <- 5

a <- as.integer(commandArgs(trailingOnly = TRUE))
chunk <- a[1]; nchunk <- a[2]

mine <- eps_all[seq(chunk, length(eps_all), by = nchunk)]
cat(sprintf("chunk %d/%d : %d endpoints\n", chunk, nchunk, length(mine)))

clip_p <- function(p) { p[p <= 0] <- 1e-300; p[p >= 1] <- 1 - 1e-16; p }

res <- vector("list", length(mine))
for (i in seq_along(mine)) {
  ep <- mine[i]
  s <- sexmap[[ep]]; if (is.null(s) || is.na(s)) s <- "both"
  tis <- TISSUES
  if (s == "female") tis <- setdiff(tis, MALE_ONLY)
  if (s == "male")   tis <- setdiff(tis, FEMALE_ONLY)

  parts <- lapply(tis, function(T) {
    f <- sprintf("%s/%s/fusios_%s_%s.txt", TS, T, ep, T)
    if (!file.exists(f)) return(NULL)
    d <- tryCatch(fread(f, select = c("gene","niv","b_fusio","se_fusio","p_fusio")),
                  error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    d <- d[niv >= MIN_IV & is.finite(b_fusio) & is.finite(se_fusio) &
           se_fusio > 0 & is.finite(p_fusio)]
    if (!nrow(d)) return(NULL)
    d[, tissue := T][]
  })
  x <- rbindlist(Filter(Negate(is.null), parts))
  if (!nrow(x)) next
  x[, gene := sub("[.].*", "", gene)]
  x[, p_clip := clip_p(p_fusio)]
  x[, tanv := tan((0.5 - p_clip) * pi)]
  x[, w := 1 / (se_fusio^2)]
  setorder(x, gene, p_fusio)
  top <- x[, .SD[1], by = gene][, .(gene, top_b = b_fusio, top_se = se_fusio,
                                    top_tissue = tissue, min_p = p_fusio)]
  agg <- x[, .(n_tissue = .N, tsum = sum(tanv), sw = sum(w), swb = sum(w * b_fusio)), by = gene]
  g <- merge(agg, top, by = "gene")
  g[, p_acat := 0.5 - atan(tsum / n_tissue) / pi]
  g[, top_z := top_b / top_se]
  g[, meta_z := (swb / sw) / sqrt(1 / sw)]
  g[, `:=`(phenocode = ep, sex_used = s, n_tis_avail = length(tis))]
  res[[i]] <- g[, .(phenocode, gene, n_tissue, n_tis_avail, sex_used,
                    p_acat, top_z, meta_z, top_tissue, min_p)]
  if (i %% 25 == 0) cat(sprintf("  [%d/%d] %s\n", i, length(mine), ep))
}
R <- rbindlist(Filter(Negate(is.null), res))
fwrite(R, sprintf("%s/acat50_chunk%03d.tsv.gz", OUT, chunk), sep = "\t")
cat(sprintf("DONE chunk %d : %d rows, %d endpoints\n", chunk, nrow(R), uniqueN(R$phenocode)))
