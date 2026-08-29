#!/usr/bin/env Rscript
## Supplementary: Figure-4B-style prevalence~coupling panels for the two classes
## the main figure does not show. A: male reproductive tract (n=39) -- the
## mirror of the female-tract pattern, direction only (Steiger P=0.55).
## B: sex-shared (n=2,492) -- near-identical slopes; the residual male-ward
## offset (Z=2.05, P=0.04) is an order of magnitude smaller than the
## gynaecological contrast and matches the sex-shared delta-r_g baseline.
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

## Fig 4B house style: coupling on x, prevalence on y (log), M solid / F dashed
panel <- function(cls, title, pt = 1.1, al = 0.55) {
  x <- s[class == cls & is.finite(rho_father) & is.finite(rho_mother)]
  dd <- rbind(data.table(lp = log10(x$prev), rho = x$rho_father, axis = "Male"),
              data.table(lp = log10(x$prev), rho = x$rho_mother, axis = "Female"))
  st <- dd[, .(r = cor(rho, lp), p = cor.test(rho, lp)$p.value), by = axis]
  sg <- steiger(x$rho_father, x$rho_mother, log10(x$prev))
  st[, lab := sprintf("'%s: r=%.2f, '*italic(P)==%s", substr(axis, 1, 1), r, pmP(p))]
  ggplot(dd, aes(rho, lp, colour = axis)) +
    geom_vline(xintercept = 0, colour = "grey85", linetype = 2, linewidth = 0.3) +
    geom_point(size = pt, alpha = al, stroke = 0) +
    geom_smooth(aes(linetype = axis), method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.8) +
    geom_text(data = st[axis == "Male"], aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.05, vjust = 1.5, size = 2.4, parse = TRUE, show.legend = FALSE) +
    geom_text(data = st[axis == "Female"], aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.05, vjust = 3, size = 2.4, parse = TRUE, show.legend = FALSE) +
    ## the Steiger contrast is reported in the figure legend, not on the panel
    scale_colour_manual(values = COLS, name = NULL,
                        labels = c(Male = "Male fertility", Female = "Female fertility")) +
    scale_linetype_manual(values = c(Male = "solid", Female = "22"), guide = "none") +
    scale_y_continuous(breaks = log10(c(0.01, 0.1, 1, 10)),
                       labels = c("0.01", "0.1", "1", "10"),
                       name = "Population prevalence (%)") +
    scale_x_continuous(name = expression("disease-fertility "*rho[GE]*" (expression effect)")) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(subtitle = sprintf("%s  (n = %s diseases)", title,
                            format(nrow(x), big.mark = ","))) +
    theme_classic(base_size = 10) +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8, colour = "#1F2225"),
          legend.position = c(0.84, 0.10), legend.text = element_text(size = 7.5),
          legend.key.size = unit(3.5, "mm"), legend.background = element_blank())
}

## the male reproductive tract panel now sits in Fig 5 as the mirror of the
## gynaecological panel, at the same 252 x 226 pt as the other Fig 5 panels
ggsave(file.path(OUT, "Fig5c_maletract.pdf"), panel("male_limited", "Male reproductive tract"),
       width = 252/72, height = 226/72, device = cairo_pdf, bg = "white")

## the supplement keeps the sex-shared class, which is the background the two
## sex-limited classes are read against
## both classes not shown in Fig 5b: the male reproductive tract and the sex-shared background
fig <- plot_grid(panel("male_limited", "Male reproductive tract"),
                 panel("sex_shared", "Sex-shared", pt = 0.55, al = 0.3),
                 ncol = 2, labels = c("a", "b"), label_size = 11, label_fontface = "bold")
ggsave(file.path(OUT, "FigS_maletract_prevalence.pdf"), fig, width = 183, height = 82,
       units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "FigS_maletract_prevalence.png"), fig, width = 183, height = 82,
       units = "mm", dpi = 300, bg = "white")
cat("written FigS_maletract_prevalence (A male tract, B sex-shared, Fig4B style)\n")
