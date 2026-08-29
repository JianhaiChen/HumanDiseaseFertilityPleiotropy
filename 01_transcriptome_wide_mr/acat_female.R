suppressPackageStartupMessages({ library(data.table) })

acat <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(NA_real_)
  p <- pmin(pmax(p, 1e-300), 1 - 1e-15)
  T_stat <- mean(tan((0.5 - p) * pi))
  if (!is.finite(T_stat) || T_stat > 1e15) return(.Machine$double.xmin)
  pa <- 0.5 - atan(T_stat) / pi
  pmin(pmax(pa, .Machine$double.xmin), 1 - 1e-15)
}

files <- list.files(".", pattern = "__female\\.txt\\.gz$", full.names = FALSE)
cat(sprintf("[INFO] %d female files\n", length(files)))

read_one <- function(f) {
  dat <- fread(f, showProgress = FALSE)
  dat <- dat[is.finite(p1) & is.finite(p2) & is.finite(b1) & is.finite(b2) &
             is.finite(se1) & is.finite(se2) & se1 > 0 & se2 > 0]
  dat[, raw_name := sub("__female\\.txt\\.gz$", "", f)]
  dat
}
ts <- proc.time()
all_dat <- rbindlist(lapply(files, read_one), use.names=TRUE, fill=TRUE)
cat(sprintf("[INFO] read %d cells in %.1fs\n", nrow(all_dat), (proc.time()-ts)[3]))
all_dat[, z1 := b1/se1][, z2 := b2/se2]

gd_acat <- all_dat[, .(
  n_tissue   = .N,
  p_acat_1   = acat(p1),
  p_acat_2   = acat(p2),
  b1_at_maxZ = b1[which.max(abs(z1))],
  b2_at_maxZ = b2[which.max(abs(z2))],
  max_abs_z1 = max(abs(z1)),
  max_abs_z2 = max(abs(z2))
), by = .(gene, raw_name, disease)]

gd_acat[, padj_acat_1 := p.adjust(p_acat_1, "bonferroni"), by = raw_name]
gd_acat[, padj_acat_2 := p.adjust(p_acat_2, "bonferroni"), by = raw_name]
gd_acat[, combo := paste0(ifelse(b1_at_maxZ >= 0, "+","-"),
                          ifelse(b2_at_maxZ >= 0, "+","-"))]
gd_acat[, class := ifelse(combo %in% c("++","--"), "Ant", "Syn")]

fwrite(gd_acat, "ACAT_female_gene_disease_all.csv")
cat(sprintf("[DONE] %d female (gene, disease) records\n", nrow(gd_acat)))
