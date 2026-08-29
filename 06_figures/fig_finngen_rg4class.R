#!/usr/bin/env Rscript
## Supplementary figure fig:fgrg4class — the four effect-direction classes against
## disease-fertility coupling across the 2,754 FinnGen diseases.
## The original generator was lost; this rebuild mirrors S10_rg4class_noband.R
## (the 97-disease counterpart, Supplementary Fig. 13) panel for panel, with the
## y range expanded at the top so the in-panel statistics clear the points.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
D   <- "${ROOT}"
OUT <- file.path(D, "fusio2trait/figures_final")

sciP <- function(p) { e <- floor(log10(p))
  paste0(sprintf("%.1f", p/10^e), "×10",
         chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(e))) }
pfmt <- function(p) fifelse(p < 1e-300, "<10⁻³⁰⁰",
                    fifelse(p < 1e-3, paste0("=", sciP(pmax(p, 1e-300))),
                            paste0("=", sprintf("%.2g", p))))

cl <- fread(file.path(D,"fusio2trait/sv_antisense_analysis_20260813/antag50_classes.tsv"))[mode=="nominal"]
st <- fread(file.path(D,"fig5_data/section_table.tsv"))[, .(phenocode, rg_father, rg_mother,
                                                           rho_father, rho_mother)]
m <- merge(cl, st, by="phenocode")
cols <- c(Male="#2C5985", Female="#C0392B")
lv   <- c("d+ f+","d- f+","d- f-","d+ f-")
BS <- 8
th <- theme_classic(base_size=BS) +
  theme(axis.title=element_text(size=BS-1), axis.text=element_text(size=BS-2, color="grey20"),
        legend.text=element_text(size=BS-1), panel.grid=element_blank())

mkrow <- function(xm, xf, xl) {
  d <- rbind(
    m[n_m>=20, .(x=get(xm), sexL="Male",   `d+ f+`=100*pp_m/n_m, `d- f+`=100*dmfp_m/n_m,
                 `d- f-`=100*mm_m/n_m, `d+ f-`=100*dpfm_m/n_m)],
    m[n_f>=20, .(x=get(xf), sexL="Female", `d+ f+`=100*pp_f/n_f, `d- f+`=100*dmfp_f/n_f,
                 `d- f-`=100*mm_f/n_f, `d+ f-`=100*dpfm_f/n_f)])
  lg <- melt(d, id.vars=c("x","sexL"), variable.name="combo", value.name="pct")[is.finite(x)]
  lg[, `:=`(combo = factor(combo, levels=lv),
            sexL  = factor(sexL, levels=c("Male","Female")))]
  rc <- lg[, {ct <- cor.test(pct, x); .(r=ct$estimate, p=ct$p.value)}, by=.(combo, sexL)]
  rc[, lab := sprintf("%s r=%.2f, P%s", fifelse(sexL=="Male","M","F"), r, pfmt(p))]
  print(rc[order(combo, sexL)])
  ggplot(lg, aes(x, pct, color=sexL, linetype=sexL, shape=sexL)) +
    ## 2,754 endpoints x 2 sexes: lighter than the 97-disease panel or the cloud goes solid
    geom_point(size=.35, alpha=.16, stroke=0) +
    geom_smooth(method="lm", formula=y~x, se=FALSE, linewidth=.6) +
    facet_wrap(~combo, nrow=1) +
    geom_text(data=rc[sexL=="Male"],   aes(x=-Inf, y=Inf, label=lab), color=cols[["Male"]],
              hjust=-.06, vjust=1.5, size=2.4, inherit.aes=FALSE) +
    geom_text(data=rc[sexL=="Female"], aes(x=-Inf, y=Inf, label=lab), color=cols[["Female"]],
              hjust=-.06, vjust=3.1, size=2.4, inherit.aes=FALSE) +
    scale_color_manual(values=cols, name=NULL) +
    scale_linetype_manual(values=c(Male="solid", Female="22"), name=NULL) +
    scale_shape_manual(values=c(Male=17, Female=4), name=NULL) +
    ## headroom for the two statistics lines
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    labs(x=xl, y="% pleiotropic genes") + th +
    theme(strip.text=element_text(face="bold", size=BS-1),
          strip.background=element_rect(fill="grey92", color=NA))
}
pR <- mkrow("rg_father","rg_mother", expression(italic(r)[g]*" (LDSC, genome-wide SNPs)"))
pG <- mkrow("rho_father","rho_mother", expression(rho[GE]*" (MR, expression effect)"))
fig <- pR/pG + plot_layout(guides="collect") + plot_annotation(tag_levels="a") &
  theme(legend.position="bottom", plot.tag=element_text(face="bold", size=12))
ggsave(file.path(OUT,"FigS_finngen_rg4class.pdf"), fig, width=180, height=118, units="mm",
       device=cairo_pdf)
ggsave(file.path(OUT,"FigS_finngen_rg4class.png"), fig, width=180, height=118, units="mm", dpi=190)
cat("written FigS_finngen_rg4class (S13 style)\n")
