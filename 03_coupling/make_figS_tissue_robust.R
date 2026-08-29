#!/usr/bin/env Rscript
## Supplementary figure fig:tissuerobust — does the antagonistic/synergistic call
## depend on taking the direction from the minimum-P tissue rather than from the
## cross-tissue IVW meta-analysis?
## Recovered 2026-08-24 from the 2026-08-14 session transcript (it had been run as an
## inline heredoc and never saved) and re-pointed at the current workbook, whose
## Supplementary Table 7 now carries the meta-analytic directions itself.
suppressMessages({library(data.table); library(ggplot2); library(openxlsx)})
D  <- "${ROOT}"
S7 <- setDT(read.xlsx(file.path(D,"fusio2trait/Supplementary Table.xlsx"),
                      sheet="Supplementary Table 7", startRow=2))
M <- S7[!is.na(pleiotropy_type_meta)]
M[, `:=`(ant_top  = pleiotropy_type      == "antagonistic (conflict)",
         ant_meta = pleiotropy_type_meta == "antagonistic (conflict)")]
cat(sprintf("pairs %d | same classification %.2f%%\n", nrow(M), 100*mean(M$ant_top==M$ant_meta)))
per <- M[, .(a=100*mean(ant_top), b=100*mean(ant_meta), n=.N),
         by=.(disease, axis=fertility_trait)][n>=10]
## The two definitions classify every pair identically, so the fractions are the
## same numbers and a correlation between them is vacuous (it is 1 by construction).
## Report the count instead.
rr <- per[, .(lab=sprintf("'%s: %d/%d diseases identical'", substr(axis,1,1),
                          sum(a==b), .N)), by=axis]
print(per[, .(diseases=.N, r=round(cor(a,b),4)), by=axis])
COLS <- c(`Male fertility`="#2C5985", `Female fertility`="#C0392B")
p <- ggplot(per, aes(a, b, colour=axis)) +
  geom_abline(intercept=0, slope=1, colour="grey80", linetype=2, linewidth=0.4) +
  geom_point(aes(size=n), alpha=0.55, stroke=0) +
  scale_size_continuous(range=c(0.8,3), name="classified pairs", breaks=c(10,30,60)) +
  scale_colour_manual(values=COLS, name=NULL) +
  geom_text(data=rr[axis=="Male fertility"],   aes(x=-Inf,y=Inf,label=lab), hjust=-0.1, vjust=1.5,
            size=2.4, parse=TRUE, show.legend=FALSE) +
  geom_text(data=rr[axis=="Female fertility"], aes(x=-Inf,y=Inf,label=lab), hjust=-0.1, vjust=3,
            size=2.4, parse=TRUE, show.legend=FALSE) +
  labs(x="antagonistic fraction, minimum-P-tissue directions (%)",
       y="antagonistic fraction,\ncross-tissue meta-analytic directions (%)") +
  theme_classic(base_size=10) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8.5),
        legend.position=c(0.99,0.02), legend.justification=c(1,0),
        legend.text=element_text(size=7), legend.key.size=unit(3,"mm"),
        legend.box="horizontal", legend.background=element_blank())
ggsave(file.path(D,"fusio2trait/figures_final/FigS_tissue_robust.pdf"), p,
       width=110, height=95, units="mm", device=cairo_pdf, bg="white")
cat("FigS_tissue_robust.pdf written\n")
