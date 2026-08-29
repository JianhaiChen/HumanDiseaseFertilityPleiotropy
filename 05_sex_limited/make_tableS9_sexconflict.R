#!/usr/bin/env Rscript
## Table S9: per-endpoint sexual-conflict summary behind Fig 4C/D and the
## Hodges-Lehmann antagonistic-share contrast.
suppressMessages(library(data.table))
D  <- "${ROOT}"
st <- fread(file.path(D, "fig5_data/section_table.tsv"))
pa <- fread(file.path(D, "fig5_data/panelD_antag.tsv"))
m <- merge(st[, .(phenocode, endpoint = pheno, class, n_cases = num_cases,
                  rho_GE_male = rho_father, rho_GE_female = rho_mother,
                  rg_male = rg_father, rg_female = rg_mother,
                  delta_rg = rg_father - rg_mother, prevalence_pct = prev,
                  prevalence_source = prev_source)],
           pa[thr == 0.05, .(phenocode, n_genes_male = n_m, n_genes_female = n_f,
                             antag_share_male_pct = pct_m,
                             antag_share_female_pct = pct_f,
                             antag_share_diff_pct = d)],
           by = "phenocode", all.x = TRUE)
## verify the manuscript numbers reproduce from this table
g  <- m[class == "gynaecological" & !is.na(delta_rg)]
ss <- m[class == "sex_shared"     & !is.na(delta_rg)]
cat(sprintf("delta_rg medians: gynae %.3f (n=%d), sex-shared %.3f (n=%d); Wilcoxon P=%.2g\n",
    median(g$delta_rg), nrow(g), median(ss$delta_rg), nrow(ss),
    wilcox.test(g$delta_rg, ss$delta_rg)$p.value))
gd <- m[class == "gynaecological", antag_share_diff_pct]; gd <- gd[!is.na(gd)]
sd_ <- m[class == "sex_shared",    antag_share_diff_pct]; sd_ <- sd_[!is.na(sd_)]
wt <- wilcox.test(gd, sd_, conf.int = TRUE)
cat(sprintf("antag share diff: HL shift %.2f pp (n=%d gynae), P=%.3f; sex-shared median %.2f pp (vs 0: P=%.2f)\n",
    wt$estimate, length(gd), wt$p.value, median(sd_), wilcox.test(sd_)$p.value))
setorder(m, class, phenocode)
num <- c("rho_GE_male","rho_GE_female","rg_male","rg_female","delta_rg",
         "antag_share_male_pct","antag_share_female_pct","antag_share_diff_pct")
m[, (num) := lapply(.SD, function(x) signif(x, 4)), .SDcols = num]
fwrite(m, file.path(D, "fusio2trait/tableS9_sexconflict_endpoints.tsv"), sep = "\t")
cat(sprintf("wrote tableS9: %d endpoints, %d columns\n", nrow(m), ncol(m)))
