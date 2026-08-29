#!/usr/bin/env Rscript
## Re-assemble Figures 4 and 5.
## Figure 4 (2026-08-13 revision):
##   A antagonism scales with coupling (fig4_panels.rds, unchanged)
##   B cross-sex fertility volcano (moved here from Fig 1C; built from
##     fig5_data/b2_meta.csv, same definition as fig1_panelC_volcano.R)
##   C female-tract coupling predicts prevalence (was B)
##   D within-endpoint paired delta r_g by endpoint class (was FigS_deltarg,
##     back in the main figure per user 2026-08-13)
## Figure 5 (recreated 2026-08-13): A r_g matrix, B r_g bins — the old Fig 4 C/D
## (fig4_panels_rg.rds), written to Fig5_2panel.pdf.
## The old panel A (four-class proportions) stays in the supplement
## (FigS_fourclass.pdf, written at the end of this script).
suppressMessages({library(data.table); library(ggplot2); library(cowplot)})
D  <- "${ROOT}"
P4 <- readRDS(file.path(D, "fusios_pipeline", "fig4_panels.rds"))
## panel A legend: title the Male/Female key "diseases" so the points read as
## per-disease strata (colour and linetype must share the name to stay merged)
suppressMessages(P4$antag <- P4$antag +
  scale_color_manual(values = c(Male = "#2C5985", Female = "#C0392B"),
                     name = "diseases") +
  scale_linetype_manual(values = c(Male = "solid", Female = "22"),
                        name = "diseases"))
## Figure 5 groups disease pairs by r_g, not rho_GE: rho_GE and the gene-effect
## signs come from the same MR estimates, so classifying the pairs with r_g -- read
## from GWAS summary statistics alone -- keeps the prediction and its test separate.
PR <- readRDS(file.path(D, "fusios_pipeline", "fig4_panels_rg.rds"))
OUT <- file.path(D, "fusio2trait", "figures_final")

## house style: P values as m x 10^e (never e-notation), italic P
pmP <- function(p) ifelse(p < 1e-3,
  sprintf("%.1f%%*%%10^%d", p / 10^floor(log10(p)), floor(log10(p))),
  sprintf("'%.2g'", p))

## ---- Figure 4 panel B: female-tract coupling vs prevalence ---------------
st2 <- fread(file.path(D, "fig5_data/section_table.tsv"))
gy <- st2[class == "gynaecological" & !is.na(rho_father) & !is.na(rho_mother) &
          !is.na(prev) & prev > 0]
gb <- melt(gy, id.vars = c("phenocode", "pheno", "prev"),
           measure.vars = c("rho_father", "rho_mother"),
           variable.name = "axis", value.name = "rho")
gb[, axis := ifelse(axis == "rho_father", "Male", "Female")]
gb[, lp := log10(prev)]
sB <- gb[, .(r = cor(rho, lp), p = cor.test(rho, lp)$p.value, n = .N), by = axis]
cat("Fig4B (female-tract, n per axis):\n"); print(sB)
COLS <- c(Male = "#2C5985", Female = "#C0392B")
sB[, lab := sprintf("'%s: r=%.2f, '*italic(P)==%s", substr(axis, 1, 1), r, pmP(p))]

## common gynaecological diseases named in the text; labelled on the paternal
## (male-fertility) points, the axis that carries the claim
pick <- data.table(
  pheno = c("Endometriosis", "Polycystic ovarian syndrome, consortium definition",
            "Ovarian cyst", "Leiomyoma of uterus",
            "Excessive, frequent and irregular menstruation",
            "Other menopausal disorders"),
  lab   = c("Endometriosis", "PCOS", "Ovarian cyst", "Uterine fibroids",
            "Irregular menstruation", "Menopausal disorders"))
exB <- merge(pick, gb[axis == "Male"], by = "pheno")

pB <- ggplot(gb, aes(rho, lp, colour = axis)) +
  geom_vline(xintercept = 0, colour = "grey85", linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.1, alpha = 0.55, stroke = 0) +
  geom_smooth(aes(linetype = axis), method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.8) +
  ## labels only, no highlight point (per user): segments point at the plain
  ## blue paternal-axis dots
  ggrepel::geom_text_repel(
    data = exB, aes(rho, lp, label = lab), size = 2.2, colour = "#1F2225",
    min.segment.length = 0, segment.size = 0.25,
    box.padding = 0.3, point.padding = 0.12, seed = 7, max.overlaps = 30,
    show.legend = FALSE) +
  geom_text(data = sB[axis == "Male"], aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 1.5, size = 2.4, parse = TRUE, show.legend = FALSE) +
  geom_text(data = sB[axis == "Female"], aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 3, size = 2.4, parse = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = COLS, name = NULL,
                      labels = c(Male = "Male fertility", Female = "Female fertility")) +
  scale_linetype_manual(values = c(Male = "solid", Female = "22"), guide = "none") +
  scale_y_continuous(breaks = log10(c(0.01, 0.1, 1, 10)),
                     labels = c("0.01", "0.1", "1", "10"),
                     name = "Population prevalence (%)") +
  scale_x_continuous(name = expression("disease-fertility "*rho[GE]*" (expression effect)")) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  ## same subtitle style as the male reproductive-tract panel it now sits beside
  labs(subtitle = sprintf("Female reproductive tract  (n = %d diseases)", nrow(gb)/2)) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        plot.subtitle = element_text(size = 8, colour = "#1F2225"),
        legend.position = c(0.84, 0.12), legend.text = element_text(size = 7.5),
        legend.key.size = unit(3.5, "mm"), legend.background = element_blank())

## ---- Figure 4 panel C: within-endpoint paired delta r_g ------------------
st <- fread(file.path(D, "fig5_data/section_table.tsv"))
st <- st[!is.na(rg_father) & !is.na(rg_mother)]
st[, del := rg_father - rg_mother]
cls_lab <- c(gynaecological = "Female\nreproductive\ntract",
             sex_shared = "Sex-shared",
             male_limited = "Male\nreproductive\ntract",
             fertility_bound = "Pregnancy &\ninfertility")
st[, cls := factor(cls_lab[class], levels = cls_lab)]
med <- st[, .(m = median(del), n = .N), by = cls]
pv <- sapply(c("gynaecological", "male_limited", "fertility_bound"), function(k)
  wilcox.test(st[class == k, del], st[class == "sex_shared", del])$p.value)
## the reference class itself: one-sample test of the shared baseline against 0
pv["sex_shared"] <- wilcox.test(st[class == "sex_shared", del])$p.value
cat("Fig4C medians/n:\n"); print(med)
cat("Fig4C wilcox vs shared (shared: vs 0):", signif(pv, 2), "\n")

## example diseases: same six gynaecological endpoints labelled in panel B
ex <- merge(pick, st, by = "pheno")
setorder(ex, cls, -del)
pann <- data.table(cls = factor(cls_lab[names(pv)], levels = levels(st$cls)),
                   p = sprintf("italic(P)==%s", pmP(pv)))
cls_cols <- setNames(c("#B4432B", "#9A9A9A", "#3B6EA5", "#C88A2E"), cls_lab)

st[, clsy := factor(cls_lab[class], levels = rev(cls_lab))]   # gynae on top
med[, clsy := factor(as.character(cls), levels = levels(st$clsy))]
ex[, clsy := factor(as.character(cls), levels = levels(st$clsy))]
pann[, clsy := factor(as.character(cls), levels = levels(st$clsy))]
## two tiers above and two below the row, walking down the delta order,
## so six labels never share a shelf with an x-neighbour
ex[, tier := rep_len(c(0.42, -0.42, 0.7, -0.7), .N), by = clsy]

set.seed(1)
pC <- ggplot(st, aes(del, clsy)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55", linewidth = 0.35) +
  ## Fig 5B house style: pale violin, tinted jitter, white box, white diamond median
  geom_violin(aes(fill = clsy), colour = NA, alpha = 0.22, scale = "width",
              width = 0.85) +
  geom_jitter(aes(colour = clsy), height = 0.16, size = 0.9, alpha = 0.2,
              stroke = 0) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.34,
               fill = "white", colour = "grey30") +
  geom_point(data = med, aes(m, clsy), shape = 23, size = 1.8, fill = "white",
             colour = "black", stroke = 0.6) +
  scale_fill_manual(values = setNames(cls_cols, cls_lab), guide = "none") +
  scale_colour_manual(values = setNames(cls_cols, cls_lab), guide = "none") +
  ## no gynae example highlights here (2026-08-13): panel C already labels
  ## the same six diseases
  geom_text(data = pann, aes(-0.315, clsy, label = p), parse = TRUE, size = 2.35,
            colour = "grey20", hjust = 0, nudge_y = -0.36) +
  scale_x_continuous(
    name = expression(paste("within-endpoint ", Delta, italic(r)[g],
                            "  (male ", "\u2212", " female fertility)"))) +
  labs(y = NULL) +
  ## direction annotation lives in panel C only; C and D share the convention
  ## x clipped to the informative range; the sex-shared tails run off the panel
  coord_cartesian(xlim = c(-0.32, 0.32), ylim = c(0.45, 4.9), clip = "on") +
  theme_classic(base_size = 10) +
  ## this panel carries the shared class axis; the forest panel beside it has none
  theme(axis.text.y = element_text(size = 7.2, colour = "black", lineheight = 0.85),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 7.5),
        axis.title.x = element_text(size = 8.5),
        plot.margin = margin(4, 2, 2, 2))

## ---- Figure 4 panel B: cross-sex fertility volcano (moved from Fig 1C) ----
## definition identical to fig1_panelC_volcano.R: MHC-excluded, either-sex
## Bonferroni on gene-level ACAT P, sign from the IVW meta effect
bv <- fread(file.path(D, "fig5_data/b2_meta.csv"))
bv[, gene := sub("[.].*$", "", gene)]
cov <- fread(file.path(D, "fig5_data/gene_coords_v19.tsv"), header = FALSE,
             col.names = c("gene", "chr", "start"))
cov[, gene := sub("[.].*$", "", gene)]
bv <- bv[!gene %in% cov[chr == "6" & start >= 25e6 & start <= 34e6, gene]]
wv <- dcast(bv, gene ~ fitness, value.var = c("meta_beta_ivw", "meta_se_ivw", "p_acat"))
wv <- wv[is.finite(meta_beta_ivw_male) & is.finite(meta_beta_ivw_female)]
wv[, dz  := meta_beta_ivw_male / meta_se_ivw_male -
            meta_beta_ivw_female / meta_se_ivw_female]
## y = the WEAKER sex's significance (2026-08-13): a gene's cross-sex signal
## is only as strong as its weaker side, so both-sex-significant genes rise
## to the top and the axis reads plainly
wv[, mlp := -log10(pmax(p_acat_male, p_acat_female))]
cat(sprintf("volcano: %d genes with -log10(weaker-sex P) >= 20 excluded from panel\n",
    wv[mlp >= 20, .N]))
## both-sex significance at FDR<0.01 (2026-08-13, per user): the conflict claim
## rests on genes significant for BOTH fertility traits — 133 genes, 39
## opposite (29.3%, vs 25.0% under both-sex Bonferroni: threshold-stable)
wv[, sig := p.adjust(p_acat_male, "BH") < 0.01 &
            p.adjust(p_acat_female, "BH") < 0.01]
wv[, grp := fifelse(!sig, "ns",
        fifelse(sign(meta_beta_ivw_male) != sign(meta_beta_ivw_female), "opp", "same"))]
nv <- wv[, .N, keyby = grp]
stopifnot(nv[grp == "opp", N] == 39, nv[grp == "same", N] == 94)
plv <- wv[mlp < 20]
VCOLS <- c(ns = "grey78", opp = "#C0392B", same = "#2C5985")
## point style matched to panels A/C: semi-transparent, no stroke
main_v <- ggplot(plv, aes(dz, mlp, colour = grp)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey40") +
  geom_point(data = plv[grp == "ns"],  size = 0.25, alpha = 0.4, stroke = 0) +
  geom_point(data = plv[grp != "ns"], size = 1.1, alpha = 0.75, stroke = 0) +
  scale_colour_manual(values = VCOLS, name = NULL,
      breaks = c("ns", "opp", "same"),
      labels = c(ns   = sprintf("Not significant in both (n=%s)", format(nv[grp=="ns",N], big.mark=",")),
                 opp  = sprintf("Opposite sign (n=%d)", nv[grp=="opp",N]),
                 same = sprintf("Same sign (n=%d)", nv[grp=="same",N])),
      guide = guide_legend(nrow = 2, byrow = TRUE,
                           override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = expression(Z[male~fertility] - Z[female~fertility]),
       y = expression(-log[10]*italic(P))) +   # defined in the caption: weaker sex
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 6.8),
        legend.key.spacing.x = unit(2, "pt"), legend.key.size = unit(8, "pt"),
        legend.margin = margin(0, 0, -4, 0),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        plot.margin = margin(4, 4, 2, 2))
barv <- data.table(grp = factor(c("opp", "same"), levels = c("opp", "same")),
                   N = c(nv[grp=="opp",N], nv[grp=="same",N]))
barv[, pct := 100 * N / sum(N)]
insetv <- ggplot(barv, aes(grp, pct, fill = grp)) +
  geom_col(width = 0.62, show.legend = FALSE, alpha = 0.85) +
  geom_text(aes(label = sprintf("%d\n%.1f%%", N, pct)), vjust = -0.25, size = 2.2,
            lineheight = 0.95) +
  scale_fill_manual(values = VCOLS[c("opp", "same")]) +
  scale_x_discrete(labels = c(opp = "Opposite\nsign", same = "Same\nsign")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_classic(base_size = 8) +
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 5.8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.background = element_blank(), panel.background = element_blank())
pV <- ggdraw(main_v) + draw_plot(insetv, x = 0.74, y = 0.50, width = 0.25, height = 0.40)

## ---- Figure 4 panels saved standalone, label-free (2026-08-13: the 6-panel
## figure is assembled in LaTeX, fusios_pipeline/fig4_6panel.tex) -------------
ggsave(file.path(OUT, "Fig4_panel_volcano.pdf"), pV, width = 89, height = 80,
       units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "Fig4_panel_gynae.pdf"), pB, width = 89, height = 80,
       units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "Fig4_panel_deltarg.pdf"), pC +
         scale_x_continuous(name = expression(Delta*italic(r)[g]*" (male "*"\u2212"*" female)")),
       width = 89, height = 80, units = "mm", device = cairo_pdf, bg = "white")

## demoted panel: antagonistic fraction ~ rho_GE -> supplement
ggsave(file.path(OUT, "FigS_antag_rho.pdf"), P4$antag, width = 120, height = 90,
       units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "FigS_antag_rho.png"), P4$antag, width = 120, height = 90,
       units = "mm", dpi = 300, bg = "white")

## ---- Figure 5: the r_g matrix + bins (old Fig 4 C/D) ----------------------
fig5 <- ggdraw(plot_grid(PR$matrix, PR$bins, ncol = 2, rel_widths = c(1.15, 1))) +
  draw_plot_label(c("A", "B"), x = c(0, 1.15 / 2.15), y = c(1, 1),
                  size = 12, fontface = "bold")
ggsave(file.path(OUT, "Fig5_2panel.pdf"), fig5, width = 183, height = 88,
       units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "Fig5_2panel.png"), fig5, width = 183, height = 88,
       units = "mm", dpi = 300, bg = "white")

## ---- old panel A (four-class proportions) -> supplement ------------------
ggsave(file.path(OUT, "FigS_fourclass.pdf"), P4$fourclass, width = 120,
       height = 90, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "FigS_fourclass.png"), P4$fourclass, width = 120,
       height = 90, units = "mm", dpi = 300, bg = "white")

cat("Fig4 = antag~rho | volcano | prev~rho (gynae) | delta rg; Fig5 = rg matrix | rg bins\n")
cat("FigS_fourclass + FigS_deltarg_antag written\n")
