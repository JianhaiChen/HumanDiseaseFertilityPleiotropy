#!/usr/bin/env Rscript
## Figure 1 panel C: cross-sex antagonism in gene-level fertility effects.
## One point per gene with FusioMR estimates for both fertility traits
## (fig5_data/b2_meta.csv). x = z(male) - z(female) on the IVW meta effect,
## y = -log10(min of the two gene-level ACAT P values). Significant = global
## Bonferroni-adjusted P < 0.05 in either sex; colour = sign agreement of the
## meta effects. Inset: opposite- vs same-sign shares among significant genes.
suppressMessages({library(data.table); library(ggplot2); library(cowplot)})
D <- "${ROOT}"
b <- fread(file.path(D, "fig5_data/b2_meta.csv"))
b[, gene := sub("[.].*$", "", gene)]
## extended MHC (chr6:25-34 Mb) excluded, as everywhere in the paper
co <- fread(file.path(D, "fig5_data/gene_coords_v19.tsv"),
            header = FALSE, col.names = c("gene", "chr", "start"))
co[, gene := sub("[.].*$", "", gene)]
mhc <- co[chr == "6" & start >= 25e6 & start <= 34e6, gene]
b <- b[!gene %in% mhc]
w <- dcast(b, gene ~ fitness, value.var = c("meta_beta_ivw", "meta_se_ivw", "p_acat"))
w <- w[is.finite(meta_beta_ivw_male) & is.finite(meta_beta_ivw_female)]
w[, zm := meta_beta_ivw_male / meta_se_ivw_male]
w[, zf := meta_beta_ivw_female / meta_se_ivw_female]
w[, dz := zm - zf]
w[, mlp := -log10(pmin(p_acat_male, p_acat_female))]
w[, sig := p.adjust(p_acat_male, "bonferroni") < 0.05 |
           p.adjust(p_acat_female, "bonferroni") < 0.05]
w[, opp := sign(meta_beta_ivw_male) != sign(meta_beta_ivw_female)]
w[, grp := fifelse(!sig, "ns", fifelse(opp, "opp", "same"))]
n <- w[, .N, keyby = grp]
stopifnot(n[grp == "opp", N] == 148, n[grp == "same", N] == 236)
cat(sprintf("ns=%d opp=%d same=%d  excluded(-log10P>=20)=%d\n",
    n[grp=="ns",N], n[grp=="opp",N], n[grp=="same",N], w[mlp >= 20, .N]))

pl <- w[mlp < 20]                    # ultra-significant genes excluded from panel
labs <- c(ns   = sprintf("Not significant (n=%s)", format(n[grp=="ns",N], big.mark=",")),
          opp  = sprintf("Significant, opposite sign (n=%d)", n[grp=="opp",N]),
          same = sprintf("Significant, same sign (n=%d)", n[grp=="same",N]))
COLS <- c(ns = "grey78", opp = "#C0392B", same = "#2C5985")

main <- ggplot(pl, aes(dz, mlp, colour = grp)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey40") +
  geom_point(data = pl[grp == "ns"],  size = 0.35, alpha = 0.45) +
  geom_point(data = pl[grp != "ns"], size = 1.1) +
  scale_colour_manual(values = COLS, labels = labs, name = NULL,
                      breaks = c("ns", "opp", "same"),
                      guide = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = expression(Z[male~fertility] - Z[female~fertility]),
       y = expression(-log[10]*"(min male/female fertility "*italic(P)*")")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 7.6),
        legend.key.spacing.x = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        axis.title = element_text(size = 10))

bar <- data.table(grp = factor(c("opp", "same"), levels = c("opp", "same")),
                  N = c(n[grp=="opp",N], n[grp=="same",N]))
bar[, pct := 100 * N / sum(N)]
inset <- ggplot(bar, aes(grp, pct, fill = grp)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n%.1f%%", N, pct)), vjust = -0.25, size = 2.5, lineheight = 0.95) +
  scale_fill_manual(values = COLS[c("opp", "same")]) +
  scale_x_discrete(labels = c(opp = "Opposite\nsign", same = "Same\nsign")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 85),
                     breaks = c(0, 25, 50, 75)) +
  ggtitle("Significant genes") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 7.5), plot.background = element_rect(colour = "grey70", linewidth = 0.3))

p <- ggdraw(main) + draw_plot(inset, x = 0.63, y = 0.52, width = 0.34, height = 0.40)
ggsave(file.path(D, "fusio2trait/figures_final/Fig1_panelC_volcano.pdf"), p,
       width = 6.4, height = 4.4)
ggsave(file.path(D, "fusio2trait/figures_final/Fig1_panelC_volcano.png"), p,
       width = 6.4, height = 4.4, dpi = 200)

## print-size variant for the composed Figure 1 (placed ~3.3 in wide in the old
## panel-A slot, so fonts here are the final print sizes)
## softer palette tuned to the figure's category colours (F47C7C / 7DA0F9)
COLS2 <- c(ns = "grey80", opp = "#E06C63", same = "#6C8FD8")
main_s <- main +
  geom_point(data = pl[grp == "ns"],  size = 0.18, alpha = 0.45) +
  geom_point(data = pl[grp != "ns"], size = 0.7) +
  scale_colour_manual(values = COLS2, name = NULL,
                      breaks = c("ns", "opp", "same"),
                      labels = c(ns   = sprintf("Not significant (n=%s)",
                                                format(n[grp=="ns",N], big.mark=",")),
                                 opp  = sprintf("Opposite sign (n=%d)", n[grp=="opp",N]),
                                 same = sprintf("Same sign (n=%d)", n[grp=="same",N])),
                      guide = guide_legend(nrow = 2, byrow = TRUE,
                                           override.aes = list(size = 1.6, alpha = 1))) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 6.0),
        legend.key.spacing.x = unit(2, "pt"), legend.key.size = unit(7, "pt"),
        legend.key.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, -4, 0),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        plot.margin = margin(4, 2, 2, 2))
## share barplot as a TRANSPARENT inset (no background box) so the volcano
## keeps the full panel width and nothing is hidden behind a white rectangle
bars_s <- inset + ggtitle(NULL) +
  scale_fill_manual(values = COLS2[c("opp", "same")]) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_classic(base_size = 7) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 5.2),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.background = element_blank(), panel.background = element_blank())
bars_s$layers[[1]]$aes_params$alpha <- 0.85
bars_s$layers[[2]]$aes_params$size <- 2
p_s <- ggdraw(main_s) + draw_plot(bars_s, x = 0.78, y = 0.15, width = 0.22, height = 0.40)
ggsave(file.path(D, "fusio2trait/figures_final/Fig1_panelC_volcano_small.pdf"), p_s,
       width = 3.6, height = 2.4, bg = "transparent")
cat("wrote Fig1_panelC_volcano.pdf/.png + _small.pdf\n")
