# =============================================================================
# Supp. Fig — Disease–fertility genetic correlation is concordant between sexes
# =============================================================================
# Male (father) vs female (mother) disease–fertility correlation across 97
# diseases, at two levels:  A = genome-wide SNP r_g (LDSC);  B = expression rho_GE.
# Both come from the finalized Fig-3 source table (fig3_disease_correlations.csv,
# written by fig3_genetic_correlation.R; 97 diseases, MHC excluded on both sides).
# Category map = Supplementary Table 1 (12-category), same as Fig 3C.
# -----------------------------------------------------------------------------
suppressMessages({library(data.table); library(ggplot2); library(readxl); library(patchwork)
  has_repel <- requireNamespace("ggrepel", quietly = TRUE)})

## ---- paths (edit DATA to your local copy of the trait-annotation tables) ----
DATA <- "/Volumes/X10Pro/data/synergistic/mr_res"     # holds Supplementary_Table_1.xlsx + _corrected_v2.csv
SRC  <- "fig3_disease_correlations.csv"                # committed alongside this script (Fig-3 source)
OUT  <- "."                                            # where to write the figure

nk <- function(x) gsub("[^a-z0-9]", "", tolower(x))

## ---- both measures from the finalized Fig-3 source (97 dz, MHC-excluded) ----
d <- fread(SRC)   # cols: dk, rhoGE_female, rhoGE_male, rgM, seM, rgF, seF

## ---- category map (identical to Fig 3C) ----
st1 <- as.data.table(read_excel(file.path(DATA, "Supplementary_Table_1.xlsx"), sheet = 1, skip = 2))
setnames(st1, 2:4, c("Trait","Type","Category"))
st1 <- st1[grepl("Disease", Type, ignore.case = TRUE), .(k = nk(Trait), Category)]
brc <- fread(file.path(DATA, "Supplementary_Table_1_corrected_v2.csv"),
             select = c("b1_disease_id","Suppl_Trait"))[, k := nk(Suppl_Trait)]
brc <- merge(brc, st1, by = "k", all.x = TRUE)
map6 <- function(c) fifelse(c == "Neuropsychiatric", "Neuropsychiatric",
  fifelse(grepl("Immune", c), "Immune",
  fifelse(c %in% c("Cardiovascular","Metabolic"), "Cardiometabolic",
  fifelse(c == "Cancers", "Cancer", fifelse(c == "Reproductive", "Reprod", "Other")))))
brc[, cat6 := map6(Category)]
dcat <- unique(brc[, .(k = nk(b1_disease_id), cat6, Suppl_Trait)])
d <- merge(d, dcat, by.x = "dk", by.y = "k", all.x = TRUE)
d[is.na(cat6), cat6 := "Other"]; d[is.na(Suppl_Trait), Suppl_Trait := dk]

cols <- c(Cardiometabolic = "#C0392B", Neuropsychiatric = "#7D3C98", Immune = "#1E8449",
          Cancer = "#E67E22", Reprod = "#C2185B", Other = "grey60")
BS <- 8
th <- theme_classic(base_size = BS) + theme(
  axis.title = element_text(size = BS), axis.text = element_text(size = BS-1, color = "grey20"),
  legend.text = element_text(size = BS-1.5), legend.title = element_blank(),
  legend.key.size = unit(3, "mm"), legend.position = "none",
  legend.background = element_rect(fill = alpha("white", .65), color = NA),
  panel.grid = element_blank())

## ---- one male-vs-female panel (leg = TRUE puts the category key bottom-right) ----
mk <- function(xc, yc, ti, leg) {
  s <- d[is.finite(get(xc)) & is.finite(get(yc))]; ct <- cor.test(s[[xc]], s[[yc]])
  conc <- 100 * mean(sign(s[[xc]]) == sign(s[[yc]])); rng <- range(c(s[[xc]], s[[yc]]))
  s[, mn := (get(xc) + get(yc)) / 2][, dv := abs(get(xc) - get(yc))]
  lab <- unique(rbind(s[order(-abs(mn))][1:6], s[order(-dv)][1:3]), by = "dk")
  p <- ggplot(s, aes(.data[[xc]], .data[[yc]])) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey85") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey85") +
    geom_abline(slope = 1, linetype = 3, color = "grey55") +
    geom_smooth(method = "lm", se = FALSE, color = "grey20", linewidth = .5) +
    geom_point(aes(color = cat6), size = 1.7, alpha = .8) + scale_color_manual(values = cols) +
    annotate("text", x = rng[1], y = rng[2], hjust = 0, vjust = 1, size = 2.8,
      label = sprintf("r = %.2f, P = %.1e\n%.0f%% same sign  (n=%d)",
                      ct$estimate, ct$p.value, conc, nrow(s))) +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = bquote("male (father) " * .(ti)), y = bquote("female (mother) " * .(ti))) + th
  if (leg) p <- p + theme(legend.position = c(.98, .02), legend.justification = c(1, 0))
  if (has_repel) p + ggrepel::geom_text_repel(data = lab, aes(label = Suppl_Trait), size = 1.9,
        color = "grey30", max.overlaps = 20, segment.size = .2, min.segment.length = 0)
  else p + geom_text(data = lab, aes(label = Suppl_Trait), size = 1.9, color = "grey30", vjust = -.6)
}
pA <- mk("rgM", "rgF", bquote(italic(r)[g] ~ "(genome-wide SNPs)"), FALSE)
pB <- mk("rhoGE_male", "rhoGE_female", bquote(rho[GE] ~ "(expression effect)"), TRUE)
fig <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
ggsave(file.path(OUT, "FigS_sex_concordance_AB.pdf"), fig, width = 180, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "FigS_sex_concordance_AB.png"), fig, width = 180, height = 95, units = "mm", dpi = 200)

for (m in list(c("rgM","rgF","r_g"), c("rhoGE_male","rhoGE_female","rho_GE"))) {
  s <- d[is.finite(get(m[1])) & is.finite(get(m[2]))]; ct <- cor.test(s[[m[1]]], s[[m[2]]])
  cat(sprintf("%-6s male-female: r=%.3f P=%.1e n=%d  %.0f%% same-sign\n", m[3], ct$estimate,
              ct$p.value, nrow(s), 100 * mean(sign(s[[m[1]]]) == sign(s[[m[2]]]))))
}
