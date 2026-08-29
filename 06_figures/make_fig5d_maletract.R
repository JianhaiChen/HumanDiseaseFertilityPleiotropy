#!/usr/bin/env Rscript
## Supplementary: Fig-5c-style prevalence~coupling panels the main figure omits.
## Fig 5d: male reproductive tract, expression-effect rho_GE -- the mirror of
## Fig 5c. Geometry and theme match Fig 5c exactly (assemble_fig4_fig5.R, pB).
## Steiger Z/P is computed and printed for the caption, not drawn on the panel.
suppressMessages({library(data.table); library(ggplot2); library(cowplot)})
D <- "${ROOT}"
OUT <- file.path(D, "fusio2trait", "figures_final")
COLS <- c(Male = "#2C5985", Female = "#C0392B")
s <- fread(file.path(D, "fig5_data/section_table.tsv"))[is.finite(prev) & prev > 0]

pmP <- function(p) ifelse(p < 1e-3,
  sprintf("%.1f%%*%%10^%d", p / 10^floor(log10(p)), floor(log10(p))),
  sprintf("'%.2g'", p))

steiger <- function(x1, x2, y) {
  n <- length(y); r13 <- cor(x1, y); r23 <- cor(x2, y); r12 <- cor(x1, x2)
  rm2 <- (r13^2 + r23^2)/2; f <- (1 - r12)/(2*(1 - rm2)); h <- (1 - f*rm2)/(1 - rm2)
  Z <- (atanh(r13) - atanh(r23)) * sqrt((n - 3)/(2*(1 - r12)*h))
  list(Z = Z, P = 2*pnorm(-abs(Z)))
}

panel <- function(cls, title, measure = "rho", pt = 1.1, al = 0.55) {
  cm <- if (measure == "rho") c("rho_father","rho_mother") else c("rg_father","rg_mother")
  x <- s[class == cls & is.finite(get(cm[1])) & is.finite(get(cm[2]))]
  dd <- rbind(data.table(lp = log10(x$prev), rho = x[[cm[1]]], axis = "Male"),
              data.table(lp = log10(x$prev), rho = x[[cm[2]]], axis = "Female"))
  st <- dd[, .(r = cor(rho, lp), p = cor.test(rho, lp)$p.value), by = axis]
  sg <- steiger(x[[cm[1]]], x[[cm[2]]], log10(x$prev))
  cat(sprintf("%-16s %-4s n=%3d  M r=%+.2f P=%.3g  F r=%+.2f P=%.3g  Steiger Z=%+.2f P=%.3g\n",
              cls, measure, nrow(x), st[axis=="Male", r], st[axis=="Male", p],
              st[axis=="Female", r], st[axis=="Female", p], sg$Z, sg$P))
  st[, lab := sprintf("'%s: r=%.2f, '*italic(P)==%s", substr(axis, 1, 1), r, pmP(p))]
  xlab <- if (measure == "rho") expression("disease-fertility "*rho[GE]*" (expression effect)") else
                                expression("disease-fertility "*r[g]*" (genome-wide)")
  ggplot(dd, aes(rho, lp, colour = axis)) +
    geom_vline(xintercept = 0, colour = "grey85", linetype = 2, linewidth = 0.3) +
    geom_point(size = pt, alpha = al, stroke = 0) +
    geom_smooth(aes(linetype = axis), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    geom_text(data = st[axis == "Male"], aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.05, vjust = 1.5, size = 2.4, parse = TRUE, show.legend = FALSE) +
    geom_text(data = st[axis == "Female"], aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.05, vjust = 3, size = 2.4, parse = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = COLS, name = NULL,
                        labels = c(Male = "Male fertility", Female = "Female fertility")) +
    scale_linetype_manual(values = c(Male = "solid", Female = "22"), guide = "none") +
    scale_y_continuous(breaks = log10(c(0.01, 0.1, 1, 10)), labels = c("0.01", "0.1", "1", "10"),
                       name = "Population prevalence (%)") +
    scale_x_continuous(name = xlab) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(subtitle = sprintf("%s  (n = %s diseases)", title, format(nrow(x), big.mark = ","))) +
    theme_classic(base_size = 10) +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8, colour = "#1F2225"),
          legend.position = c(0.84, 0.12), legend.text = element_text(size = 7.5),
          legend.key.size = unit(3.5, "mm"), legend.background = element_blank())
}

fig <- panel("male_limited", "Male reproductive tract", pt = 1.3, al = 0.6)
ggsave(file.path(OUT, "Fig5d_maletract.pdf"), fig, width = 89, height = 80,
       units = "mm", device = cairo_pdf, bg = "white")
cat("written Fig5d_maletract.pdf (male reproductive tract, rho_GE)\n")
