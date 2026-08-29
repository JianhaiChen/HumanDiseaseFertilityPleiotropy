#!/usr/bin/env Rscript
## Fig 2D: pooled ANTAGONISTIC share (d+f+ + d-f-, the definitional classes)
## among disease-fertility shared genes vs log10 prevalence (50-tissue layer).
## Supp (FigS_classprev): the four directional classes resolved individually.
suppressMessages({library(data.table); library(ggplot2)})
D <- "${ROOT}"
cl <- fread(file.path(D,"fusio2trait/sv_antisense_analysis_20260813/antag50_classes.tsv"))[mode=="nominal"]
st <- fread(file.path(D,"fig5_data/section_table.tsv"))[, .(phenocode, prev)]
st <- st[!is.na(prev) & prev > 0]
m <- merge(cl, st, by="phenocode"); m[, lp := log10(prev)]
COLS <- c(Male="#2C5985", Female="#C0392B")
pmP <- function(p) ifelse(p < 1e-3,
  sprintf("%.1f%%*%%10^%d", p/10^floor(log10(p)), floor(log10(p))), sprintf("'%.2g'", p))

## ---- main panel: pooled antagonistic share --------------------------------
long <- rbind(
  m[n_m>=20, .(phenocode, lp, axis="Male",   sh = 100*(pp_m+mm_m)/n_m)],
  m[n_f>=20, .(phenocode, lp, axis="Female", sh = 100*(pp_f+mm_f)/n_f)])
## same panel for the pooled antagonistic share and for the d+f- share (Fig 4b),
## so the two panels of Fig 4 cannot drift apart in style or in the fitted numbers.
mkpanel <- function(long, xlab, vline = NA_real_, legpos = c(0.84, 0.12)) {
sB <- long[, {ct <- cor.test(sh, lp); .(r = ct$estimate, p = ct$p.value, n = .N)}, by=axis]
print(sB)
sB[, lab := sprintf("'%s: r=%.2f, '*italic(P)==%s", substr(axis,1,1), r, pmP(p))]
ggplot(long, aes(sh, lp, colour=axis)) +
  ## the 50% reference only means something for the pooled antagonistic share
  {if (is.finite(vline)) geom_vline(xintercept = vline, colour="grey85", linetype=2, linewidth=0.3)} +
  geom_point(size=0.7, alpha=0.35, stroke=0) +
  geom_smooth(aes(linetype=axis), method="lm", formula=y~x, se=FALSE, linewidth=0.8) +
  geom_text(data=sB[axis=="Male"],   aes(x=-Inf, y=Inf, label=lab), hjust=-0.05, vjust=1.5,
            size=2.4, parse=TRUE, show.legend=FALSE) +
  geom_text(data=sB[axis=="Female"], aes(x=-Inf, y=Inf, label=lab), hjust=-0.05, vjust=3,
            size=2.4, parse=TRUE, show.legend=FALSE) +
  scale_colour_manual(values=COLS, name=NULL,
                      labels=c(Male="Male fertility", Female="Female fertility")) +
  scale_linetype_manual(values=c(Male="solid", Female="22"), guide="none") +
  scale_y_continuous(breaks=log10(c(0.01,0.1,1,10)), labels=c("0.01","0.1","1","10"),
                     name="Population prevalence (%)") +
  scale_x_continuous(name=xlab) +
  guides(colour=guide_legend(override.aes=list(size=2, alpha=1))) +
  theme_classic(base_size=10) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=9),
        legend.position=legpos, legend.text=element_text(size=7.5),
        legend.key.size=unit(3.5,"mm"), legend.background=element_blank())
}
pD <- mkpanel(long, "disease-fertility antagonistic pleiotropic genes (%)", vline = 50)
long_dpfm <- rbind(
  m[n_m>=20, .(phenocode, lp, axis="Male",   sh = 100*dpfm_m/n_m)],
  m[n_f>=20, .(phenocode, lp, axis="Female", sh = 100*dpfm_f/n_f)])
pB <- mkpanel(long_dpfm, expression(paste(italic(d)["+"]*italic(f)["-"], " synergistic genes (%)")),
              legpos = c(0.19, 0.12))
OUT <- file.path(D, "fusio2trait/figures_final")
ggsave(file.path(OUT,"Fig4b_dpfm.pdf"), pB, width=89, height=80, units="mm",
       device=cairo_pdf, bg="white")
ggsave(file.path(OUT,"Fig4b_dpfm.png"), pB, width=91, height=82, units="mm",
       dpi=300, bg="white")
ggsave(file.path(OUT,"Fig2D_antagclass.pdf"), pD, width=89, height=80, units="mm",
       device=cairo_pdf, bg="white")
ggsave(file.path(OUT,"Fig2D_antagclass.png"), pD, width=91, height=82, units="mm",
       dpi=300, bg="white")

## ---- supp: the four classes resolved --------------------------------------
lv <- c("d+f+","d-f-","d+f-","d-f+")
cls <- rbind(
  m[n_m>=20, .(phenocode, lp, axis="Male",   `d+f+`=100*pp_m/n_m, `d-f-`=100*mm_m/n_m,
               `d+f-`=100*dpfm_m/n_m, `d-f+`=100*dmfp_m/n_m)],
  m[n_f>=20, .(phenocode, lp, axis="Female", `d+f+`=100*pp_f/n_f, `d-f-`=100*mm_f/n_f,
               `d+f-`=100*dpfm_f/n_f, `d-f+`=100*dmfp_f/n_f)])
lg <- melt(cls, id.vars=c("phenocode","lp","axis"), variable.name="class", value.name="sh")
lg[, class := factor(class, levels=lv)]
sS <- lg[, {ct <- cor.test(sh, lp); .(r = ct$estimate, p = ct$p.value)}, by=.(class, axis)]
print(sS)
sS[, lab := sprintf("'%s: r=%+.2f, '*italic(P)==%s", substr(axis,1,1), r, pmP(p))]
sS[, vj := fifelse(axis=="Male", 1.4, 2.9)]
pS <- ggplot(lg, aes(sh, lp, colour=axis)) +
  geom_point(size=0.5, alpha=0.3, stroke=0) +
  geom_smooth(aes(linetype=axis), method="lm", formula=y~x, se=FALSE, linewidth=0.7) +
  geom_text(data=sS, aes(x=Inf, y=Inf, label=lab, vjust=vj), hjust=1.05,
            size=2.3, parse=TRUE, show.legend=FALSE) +
  facet_wrap(~class, nrow=2, scales="free_x") +
  scale_colour_manual(values=COLS, name=NULL,
                      labels=c(Male="Male fertility", Female="Female fertility")) +
  scale_linetype_manual(values=c(Male="solid", Female="22"), guide="none") +
  scale_y_continuous(breaks=log10(c(0.01,0.1,1,10)), labels=c("0.01","0.1","1","10"),
                     name="Population prevalence (%)") +
  scale_x_continuous(name="class share among shared genes (%)") +
  guides(colour=guide_legend(override.aes=list(size=2, alpha=1))) +
  theme_bw(base_size=10) +
  theme(panel.grid=element_blank(), strip.background=element_blank(),
        strip.text=element_text(size=9, face="bold"),
        axis.text=element_text(size=7.5), axis.title=element_text(size=9),
        legend.position="top", legend.text=element_text(size=8))
ggsave(file.path(OUT,"FigS_classprev.pdf"), pS, width=160, height=150, units="mm",
       device=cairo_pdf, bg="white")
ggsave(file.path(OUT,"FigS_classprev.png"), pS, width=160, height=150, units="mm",
       dpi=300, bg="white")
cat("written Fig2D (pooled antag) + FigS_classprev (4 classes)\n")
