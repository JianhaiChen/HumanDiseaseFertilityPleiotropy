#!/usr/bin/env Rscript
## Supplementary figure + tables 11a/11b: delta r_g by disease category.
## Sex-shared endpoints only, so no reproductive-tract disease can carry a category.
## Each category vs all remaining sex-shared endpoints; BH across the 19 categories.
suppressMessages({library(data.table); library(ggplot2)})

D   <- "${ROOT}"
OUT <- file.path(D, "fusio2trait", "figures_final")

st <- fread(file.path(D, "fig5_data/section_table.tsv"))[is.finite(rg_father) & is.finite(rg_mother)]
st[, d := rg_father - rg_mother]
ss <- st[class == "sex_shared"]

r <- ss[, {
  if (.N >= 20) .(n = .N, med = median(d), lo = NA_real_, hi = NA_real_,
                  pct = 100 * mean(d > 0),
                  P = wilcox.test(d, ss[cat != .BY$cat, d])$p.value) else NULL
}, by = cat]
## Hodges-Lehmann interval for the shift against the rest of the sex-shared set
for (k in seq_len(nrow(r))) {
  w <- wilcox.test(ss[cat == r$cat[k], d], ss[cat != r$cat[k], d], conf.int = TRUE)
  set(r, k, "lo", w$conf.int[1]); set(r, k, "hi", w$conf.int[2])
  set(r, k, "med", unname(w$estimate))            # HL shift, not raw median
}
r[, q := p.adjust(P, "BH")]

## short display names for the FinnGen category strings
short <- function(x) {
  x <- sub("^[IVXLC]+ ", "", x)
  x <- sub(" \\([A-Za-z0-9_]+\\)$", "", x)
  x <- sub(" and connective tissue", "", x)
  x <- sub(" and subcutaneous tissue", "", x)
  x <- sub(" and blood-forming organs and certain disorders involving the immune mechanism",
           " and immune mechanism", x)
  x <- sub("Psychiatric endpoints from Katri R.*", "Psychiatric endpoints", x)
  x <- sub(", deformations and chromosomal abnormalities", "", x)
  x <- sub(", from cancer register.*", " (cancer register)", x)
  x <- sub(" from hospital discharges", " (hospital)", x)
  substr(x, 1, 44)
}
r[, lab := short(cat)]
r[, sig := fifelse(q < 0.05, "FDR < 0.05", "not significant")]
r <- r[order(med)][, lab := factor(lab, levels = lab)]

base <- median(ss$d)
gyn  <- st[class == "gynaecological", d]
gynHL <- unname(wilcox.test(gyn, ss$d, conf.int = TRUE)$estimate)

p <- ggplot(r, aes(med, lab)) +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey75") +
  geom_vline(xintercept = gynHL, linetype = 2, linewidth = .4, colour = "#C0392B") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .4, colour = "grey55") +
  geom_point(aes(size = n, colour = sig)) +
  scale_colour_manual(values = c("FDR < 0.05" = "#2C5985", "not significant" = "grey65"),
                      name = NULL) +
  scale_size_continuous(range = c(1, 3.4), name = "endpoints", breaks = c(25, 75, 150)) +
  labs(x = expression("shift in "*Delta*italic(r)[g]*" toward male (right) or female (left) fertility coupling"),
       y = NULL) +
  theme_classic(base_size = 7) +
  theme(legend.position = "right", legend.key.height = unit(8, "pt"),
        axis.text.y = element_text(size = 6.2))

ggsave(file.path(OUT, "FigS_category_deltarg.pdf"), p, width = 150, height = 78,
       units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "FigS_category_deltarg.png"), p, width = 150, height = 78,
       units = "mm", dpi = 200)

## coupling level per category (same sex-shared set), for Supplementary Table 11a
lev <- ss[, .(rho_male = median(rho_father, na.rm = TRUE), rho_female = median(rho_mother, na.rm = TRUE),
              rg_male  = median(rg_father,  na.rm = TRUE), rg_female  = median(rg_mother,  na.rm = TRUE),
              pct_rho_positive = 100 * mean(rho_father > 0, na.rm = TRUE)), by = cat]
r <- merge(r, lev, by = "cat")

fwrite(r[order(-med), .(category = cat, n_endpoints = n,
                        median_rho_GE_male = round(rho_male, 5), median_rho_GE_female = round(rho_female, 5),
                        median_rg_male = round(rg_male, 4), median_rg_female = round(rg_female, 4),
                        pct_rho_positive = round(pct_rho_positive, 1),
                        HL_shift_delta_rg = round(med, 4),
                        CI_low = round(lo, 4), CI_high = round(hi, 4),
                        pct_male_leaning = round(pct, 1),
                        P = signif(P, 3), q_BH = signif(q, 3))],
       file.path(D, "fusio2trait", "tableS11_category_delta_rg.tsv"), sep = "\t")

cat(sprintf("sex-shared endpoints %d in %d categories; baseline median delta_rg %+.4f; gynaecological HL %+.4f\n",
            nrow(ss), nrow(r), base, gynHL))
print(r[order(-med), .(lab, n, HL = round(med, 4), pct = round(pct), q = signif(q, 2))])

## Per-endpoint detail for the FDR-significant categories (Supplementary Table 11b)
sig_cats <- r[q < 0.05, as.character(cat)]
det <- ss[cat %in% sig_cats][order(cat, -d)]
fwrite(det[, .(category = cat, phenocode, endpoint = pheno, n_cases = num_cases,
               prevalence_pct = round(prev, 3), prevalence_source = prev_source,
               rg_male_fertility = round(rg_father, 4), rg_female_fertility = round(rg_mother, 4),
               delta_rg = round(d, 4), male_leaning = d > 0)],
       file.path(D, "fusio2trait", "tableS11b_category_endpoints.tsv"), sep = "\t")
cat(sprintf("per-endpoint detail: %d endpoints in %d FDR-significant categories\n", nrow(det), length(sig_cats)))
