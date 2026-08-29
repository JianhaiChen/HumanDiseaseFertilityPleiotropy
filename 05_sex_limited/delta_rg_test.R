#!/usr/bin/env Rscript
# Per-endpoint test of Delta r_g = rg(male fertility) - rg(female fertility).
# Conservative: assumes cov(rg_father_hat, rg_mother_hat) = 0, which overstates
# the SE of the difference because both estimates share the disease GWAS
# (and the two fertility GWAS are positively genetically correlated).
suppressMessages(library(data.table))

rg <- fread("${ROOT}/fig5_data/fg_rg_full.tsv")
st <- fread("${ROOT}/fig5_data/section_table.tsv")
d <- merge(rg, st[, .(phenocode, pheno, class, num_cases, prev)], by = "phenocode")
d <- d[!is.na(rg_father) & !is.na(rg_mother) & se_father > 0 & se_mother > 0]
# usable precision only: LDSC rg on tiny endpoints has huge SEs
d <- d[se_father < 0.5 & se_mother < 0.5]

d[, delta := rg_father - rg_mother]
d[, se_d  := sqrt(se_father^2 + se_mother^2)]
d[, z     := delta / se_d]
d[, p     := 2 * pnorm(-abs(z))]
d[, q     := p.adjust(p, "BH")]

cat(sprintf("endpoints tested: %d\n", nrow(d)))
sig <- d[q < 0.05][order(-abs(z))]
cat(sprintf("FDR<0.05: %d  (male-leaning %d / female-leaning %d)\n\n",
            nrow(sig), sig[delta > 0, .N], sig[delta < 0, .N]))

cat("== class composition of significant endpoints ==\n")
print(sig[, .N, by = .(class, dir = ifelse(delta > 0, "male", "female"))])

cat("(conservative pass done; calibrated pass below)\n")

## ---- empirical calibration ----------------------------------------------
# cov=0 overstates se_d, so z is deflated by a common factor. Estimate the
# deflation from the sex-shared endpoints (empirical null, centred on their
# own median so the small shared baseline does not inflate it; MAD for
# robustness against true outliers).
lam <- d[class == "sex_shared", mad(z, center = median(z))]
cat(sprintf("\nempirical z scale among sex-shared (MAD): %.3f\n", lam))
d[, z_adj := (z - 0) / lam]
d[, p_adj := 2 * pnorm(-abs(z_adj))]
d[, q_adj := p.adjust(p_adj, "BH")]

sig <- d[q_adj < 0.05][order(-abs(z_adj))]
cat(sprintf("calibrated FDR<0.05: %d  (male-leaning %d / female-leaning %d)\n\n",
            nrow(sig), sig[delta > 0, .N], sig[delta < 0, .N]))
print(sig[, .N, by = .(class, dir = ifelse(delta > 0, "male", "female"))])

if (nrow(sig)) {
  cat("\n== calibrated significant MALE-leaning (top 25) ==\n")
  print(sig[delta > 0][1:min(25, sum(sig$delta > 0)),
        .(pheno = substr(pheno, 1, 50), class, delta = round(delta, 3),
          z_adj = round(z_adj, 2), q_adj = signif(q_adj, 2))])
  cat("\n== calibrated significant FEMALE-leaning (top 25) ==\n")
  print(sig[delta < 0][1:min(25, sum(sig$delta < 0)),
        .(pheno = substr(pheno, 1, 50), class, delta = round(delta, 3),
          z_adj = round(z_adj, 2), q_adj = signif(q_adj, 2))])
  tab <- table(gynae = d$class == "gynaecological", sig_male = d$q_adj < 0.05 & d$delta > 0)
  cat("\n"); print(fisher.test(tab))
}
fwrite(d[order(q_adj)], "${ROOT}/fig5_data/delta_rg_test.tsv", sep = "\t")
