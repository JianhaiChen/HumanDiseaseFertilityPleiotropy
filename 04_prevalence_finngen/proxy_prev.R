## Can FinnGen's own case fraction stand in for FinRegistry prevalence, so the 144
## cancer endpoints can enter Figure 2C at all? FinRegistry publishes no prevalence
## for cancer-registry endpoints (115/116 flagged no_summary_statistics), so the
## alternative is to leave every cancer out of a claim made about disease generally.
suppressMessages(library(data.table))
man <- fread("${BASEDIR}/finngen_r13/R13_manifest.tsv",
             select = c("phenocode", "num_cases", "num_controls"))
pv  <- fread("${BASEDIR}/prevalence_full.tsv")
man[, cf := 100 * num_cases / (num_cases + num_controls)]
d <- merge(man, pv, by = "phenocode")[cf > 0 & finregistry_prev_all > 0]

cat(sprintf("endpoints with both measures: %d\n", nrow(d)))
ct <- cor.test(log10(d$cf), log10(d$finregistry_prev_all))
cat(sprintf("  agreement on the log scale: r=%.3f, slope=%.3f\n",
            ct$estimate, coef(lm(log10(finregistry_prev_all) ~ log10(cf), d))[2]))
cat(sprintf("  median ratio (registry/cohort): %.2f\n",
            median(d$finregistry_prev_all / d$cf)))

ca <- man[grepl("^C3_", phenocode)]
cat(sprintf("\ncancer endpoints: %d, case fraction median %.3f%% (range %.4f-%.2f)\n",
            nrow(ca), median(ca$cf), min(ca$cf), max(ca$cf)))
cat(sprintf("non-cancer with registry prevalence: case fraction median %.3f%%\n",
            d[!grepl("^C3_", phenocode), median(cf)]))

## does the Figure 2C conclusion hold when the proxy replaces the registry value,
## and does it survive adding the cancers?
rho <- fread("${BASEDIR}/rhoGE_50tissue_matched.tsv")
r <- merge(rho[!is.na(rho_father), .(phenocode, rho_father, rho_mother)], man, by = "phenocode")
r <- merge(r, pv, by = "phenocode", all.x = TRUE)
r[, is_ca := grepl("^C3_", phenocode)]

cat("\n=== Figure 2C under three definitions ===\n")
f <- function(x, xv, lab) {
  x <- x[is.finite(get(xv)) & get(xv) > 0]
  if (nrow(x) < 50) { cat(sprintf("  %-40s n=%d too few\n", lab, nrow(x))); return(invisible()) }
  a <- cor.test(x$rho_father, log10(x[[xv]])); b <- cor.test(x$rho_mother, log10(x[[xv]]))
  cat(sprintf("  %-40s n=%-5d male r=%+.3f (P=%.1e) | female r=%+.3f (P=%.1e)\n",
              lab, nrow(x), a$estimate, a$p.value, b$estimate, b$p.value))
}
cat("  published                                n=1,786 male r=+0.340 | female r=+0.230\n")
f(r, "finregistry_prev_all", "registry prevalence (no cancers)")
f(r[is_ca == FALSE], "cf", "cohort case fraction, cancers excluded")
f(r, "cf", "cohort case fraction, cancers included")
f(r[is_ca == TRUE], "cf", "cancers only")
