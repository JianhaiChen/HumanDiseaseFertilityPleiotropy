#!/usr/bin/env Rscript
## Two supplementary counterparts, so each coupling measure is shown for both analyses:
##   a) delta rho_GE by disease class  -- expression-layer counterpart of Fig. 5c
##   b) prevalence by r_g sign         -- SNP-layer counterpart of the rho_GE violin
## Same conventions as Fig. 5c/d: non-reference classes vs sex-shared, sex-shared
## itself one-sample against zero.
suppressMessages({library(data.table); library(ggplot2)})

D   <- "${ROOT}"
OUT <- file.path(D, "fusio2trait", "figures_final")
d   <- fread(file.path(D, "fusio2trait", "tableS9_sexconflict_endpoints.tsv"))[class != ""]
d[, delta_rho := rho_GE_male - rho_GE_female]

LAB  <- c(gynaecological = "Female\nreproductive\ntract", sex_shared = "Sex-shared",
          male_limited = "Male\nreproductive\ntract", fertility_bound = "Pregnancy &\ninfertility")
COLS <- setNames(c("#B4432B", "#9A9A9A", "#3B6EA5", "#C88A2E"), LAB)
pmP <- function(p) {
  if (p >= 0.001) return(sprintf("italic(P)==%s", format(signif(p, 2), scientific = FALSE)))
  e <- floor(log10(p)); sprintf("italic(P)==%s%%*%%10^%d", signif(p / 10^e, 2), e)
}

x   <- d[is.finite(delta_rho), .(v = delta_rho, class)]
ref <- x[class == "sex_shared", v]
st  <- x[, .(cls = factor(LAB[class], levels = rev(LAB)), v)]
med <- st[, .(m = median(v)), by = cls]
## one-sample signed-rank against zero for every class, matching Fig. 6c/d.
pv  <- x[, .(p = wilcox.test(v)$p.value), by = class]
pv[, `:=`(cls = factor(LAB[class], levels = levels(st$cls)), lab = vapply(p, pmP, ""))]

set.seed(1)
pA <- ggplot(st, aes(v, cls)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55", linewidth = 0.35) +
  geom_violin(aes(fill = cls), colour = NA, alpha = 0.22, scale = "width", width = 0.85) +
  geom_jitter(aes(colour = cls), height = 0.15, size = 0.55, alpha = 0.30, stroke = 0) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.34, fill = "white", colour = "grey30") +
  geom_point(data = med, aes(m, cls), shape = 23, size = 1.7, fill = "white", colour = "black", stroke = 0.55) +
  geom_text(data = pv, aes(-0.049, cls, label = lab), parse = TRUE, size = 2.2,
            colour = "grey20", hjust = 0, nudge_y = -0.36) +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_colour_manual(values = COLS, guide = "none") +
  scale_x_continuous(name = expression(Delta * rho[GE] * "  (male minus female fertility)")) +
  labs(y = NULL) +
  coord_cartesian(xlim = c(-0.05, 0.05), ylim = c(0.45, 4.9), clip = "on") +
  theme_classic(base_size = 7.4) +
  theme(axis.text.y = element_text(size = 6.6, lineheight = 0.95),
        plot.margin = margin(4, 6, 3, 3))
ggsave(file.path(OUT, "FigS_deltarho_class.pdf"), pA, width = 252/72, height = 226/72, device = cairo_pdf)

## ---- prevalence by coupling sign, both measures as panels a and b ----------
suppressMessages(library(patchwork))
viol <- function(meas, xlab) {
  q <- melt(d[prevalence_pct > 0], id.vars = c("phenocode", "prevalence_pct"),
            measure.vars = paste0(meas, c("_male", "_female")),
            variable.name = "axis", value.name = "v")[is.finite(v)]
  q[, axis := factor(fifelse(grepl("_male$", axis), "Male", "Female"), levels = c("Male", "Female"))]
  q[, grp := factor(fifelse(v > 0, "pos", "neg"), levels = c("neg", "pos"))]
  sv <- q[, .(pp = wilcox.test(prevalence_pct[grp == "pos"], prevalence_pct[grp == "neg"])$p.value), by = axis]
  mq <- q[, .(m = median(prevalence_pct)), by = .(axis, grp)]
  set.seed(1)
  ggplot(q, aes(grp, prevalence_pct)) +
    geom_violin(aes(fill = grp), colour = NA, alpha = 0.24, scale = "width", width = 0.85) +
    geom_jitter(aes(colour = grp), width = 0.28, height = 0, size = 0.28, alpha = 0.16, stroke = 0) +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.35, fill = "white", colour = "grey30") +
    geom_point(data = mq, aes(grp, m), shape = 23, size = 1.5, fill = "white", colour = "black", stroke = 0.5) +
    facet_wrap(~axis) +
    geom_text(data = sv, aes(1.5, 60, label = sprintf("Mann-Whitney~%s", vapply(pp, pmP, ""))),
              parse = TRUE, size = 2.3, inherit.aes = FALSE) +
    scale_fill_manual(values = c(neg = "#2C5985", pos = "#C0392B"), guide = "none") +
    scale_colour_manual(values = c(neg = "#2C5985", pos = "#C0392B"), guide = "none") +
    scale_x_discrete(labels = c(neg = "< 0", pos = "> 0")) +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
    labs(x = xlab, y = "Population prevalence (%)") +
    theme_bw(7.4) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold", size = 6.6))
}
pB <- viol("rho_GE", expression("disease-fertility " * rho[GE])) /
      viol("rg",     expression("disease-fertility " * italic(r)[g])) +
      plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 10))
ggsave(file.path(OUT, "FigS_prevalence_violin_both.pdf"), pB, width = 130, height = 155,
       units = "mm", device = cairo_pdf, bg = "white")

print(pv[, .(class, P = signif(p, 3))])
