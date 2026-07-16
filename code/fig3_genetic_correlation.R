# =============================================================================
# Figure 3 — Disease–fertility genetic correlation and antagonistic pleiotropy
# =============================================================================
# Two complementary estimators of the per-disease disease–fertility genetic
# correlation, computed SEX-SPECIFICALLY (father = male, mother = female):
#   rho_GE : gene-expression correlation. cor(z_disease, z_fitness) across genes,
#            where z = MR effect (top_beta / top_se) from FusioMR (multi-tissue ACAT).
#   r_g    : genome-wide SNP genetic correlation from LD-score regression (LDSC).
#
# ---------------------------------------------------------------------------
# MHC HANDLING  (chr6:25–34 Mb) — both estimators EXCLUDE MHC (consistent)
# ---------------------------------------------------------------------------
#   * r_g: the standard European LD-score reference used here (eur_w_ld_chr)
#     contains NO SNPs in chr6:24,999,740–34,005,069 (a 9-Mb gap = the MHC).
#     LDSC --rg can only use SNPs present in BOTH the sumstats and the reference,
#     so MHC SNPs are automatically dropped -> our r_g ALREADY EXCLUDES MHC.
#     (Verified directly from eur_w_ld_chr/6.l2.ldscore.gz.)
#   * rho_GE: excludes MHC by the explicit !isMHC(chr,pos) filter below.
#   * Therefore BOTH estimators exclude MHC -> the comparison is on the same
#     footing. INCLUDE_MHC <- FALSE is the correct/consistent setting.
#   * Rationale (see Methods): the MHC's extreme long-range LD prevents reliable
#     causal-gene identification in gene-level MR/TWAS (LD-induced co-regulation
#     tags many non-causal genes) and biases SNP-level LD-score regression;
#     hence it is removed as a region. MHC genes are nonetheless retained and
#     flagged (in_MHC) in Supplementary Table 2 — exclusion is a causal-resolution
#     limitation, not a claim the MHC is biologically unimportant.
#   * Sensitivity: for rho_GE, including MHC (as pseudo-replicated genes) or as
#     ~7 recombination-block representatives changes per-disease rho_GE by <=0.02
#     (representative vs excluded: male cor 1.00, female 1.00) — negligible.
# =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

## ---- paths (edit DATA to your local copy) ----
DATA <- "/Volumes/X10Pro/data/synergistic/mr_res"
OUT  <- "."                                   # where to write Fig3 pdf/png
INCLUDE_MHC <- FALSE                           # see MHC HANDLING note above
Z_CAP <- 9                                    # drop |z|>9 outlier genes (robustness)

strip <- function(x) sub("[.].*", "", x)
nk    <- function(x) gsub("[^a-z0-9]", "", tolower(x))
isMHC <- function(chr, pos) chr == 6 & pos > 25e6 & pos < 34e6

co <- unique(fread(file.path(DATA, "gene_coords_v19.tsv"),
                   col.names = c("g","chr","pos"))[, g := strip(g)], by = "g")

## ---- rho_GE : per disease x sex, gene-expression correlation (MR) ----
b1 <- fread(file.path(DATA, "gene_level_b1_disease_acat_summary.dropbox.csv"),
            select = c("gene","disease","top_beta","top_se"))[
            , g := strip(gene)][is.finite(top_beta) & top_se > 0][, z1 := top_beta/top_se]
b2 <- fread(file.path(DATA, "gene_level_b2_fitness_acat_summary.dropbox.csv"),
            select = c("gene","fitness","top_beta","top_se"))[
            , g := strip(gene)][is.finite(top_beta) & top_se > 0][, z2 := top_beta/top_se]
b1 <- merge(b1, co, by = "g"); b2 <- merge(b2, co, by = "g")
keep <- function(chr, pos) if (INCLUDE_MHC) TRUE else !isMHC(chr, pos)
za <- b1[keep(chr, pos) & abs(z1) < Z_CAP]
zb <- b2[keep(chr, pos) & abs(z2) < Z_CAP]

rho <- rbindlist(lapply(unique(za$disease), function(dz)
  rbindlist(lapply(c("male","female"), function(sx) {
    m <- merge(za[disease == dz, .(g, z1)], zb[fitness == sx, .(g, z2)], by = "g")
    if (nrow(m) < 30) return(NULL)
    ct <- cor.test(m$z1, m$z2)
    data.table(dk = nk(dz), sex = sx, rge = unname(ct$estimate),
               p = ct$p.value, ng = nrow(m))
  }))))
rhoW <- dcast(rho, dk ~ sex, value.var = "rge")

## ---- r_g : LDSC (father=male, mother=female). Excludes MHC via LD-score ref. ----
Lr <- readLines(file.path(DATA, "dz_rg_97_ALLESTIMABLE.txt"))
gv <- function(l, t) {  # returns c(estimate, se)
  m <- regmatches(l, regexpr(paste0(t, "=[-0-9.eE]+ \\([0-9.eE]+\\)"), l))
  if (!length(m)) return(c(NA, NA))
  as.numeric(c(sub(paste0(t, "=([-0-9.eE]+).*"), "\\1", m),
               sub(".*\\(([0-9.eE]+)\\)", "\\1", m)))
}
rg <- rbindlist(lapply(Lr, function(l) {
  fa <- gv(l, "father"); mo <- gv(l, "mother")
  data.table(dk = nk(sub("_$", "", strsplit(l, "\t")[[1]][1])),
             rgM = fa[1], seM = fa[2], rgF = mo[1], seF = mo[2])
}))[is.finite(rgM) & is.finite(rgF)]

## ---- write the per-disease correlation table (Fig 3 source data) ----
tab <- merge(rhoW, rg, by = "dk", all = TRUE)
setnames(tab, c("male","female"), c("rhoGE_male","rhoGE_female"), skip_absent = TRUE)
fwrite(tab, file.path(OUT, "fig3_disease_correlations.csv"))

## ---- sex-specific convergence rho_GE vs r_g (both EXCLUDE MHC = same scale) ----
cv <- merge(rhoW, rg, by = "dk")
for (s in c("M","F")) {
  x <- cv[[c(M="male", F="female")[s]]]; y <- cv[[paste0("rg", s)]]
  ok <- is.finite(x) & is.finite(y); ct <- cor.test(x[ok], y[ok])
  cat(sprintf("convergence rho_GE vs r_g  %s: r=%.2f  P=%.1e  n=%d  (%.0f%% sign-concordant)\n",
              s, ct$estimate, ct$p.value, sum(ok),
              100*mean(sign(x[ok]) == sign(y[ok]))))
}
cat(sprintf("MHC included: %s  |  diseases: rho_GE %d, r_g %d\n",
            INCLUDE_MHC, nrow(rhoW), nrow(rg)))

# (Plotting of panels A/B/C is in fig3_plot.R, which reads
#  fig3_disease_correlations.csv written above.)
