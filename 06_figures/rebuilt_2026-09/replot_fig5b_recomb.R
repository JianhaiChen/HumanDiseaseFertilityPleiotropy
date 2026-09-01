#!/usr/bin/env Rscript
## Replot Fig 5b from the cached recomb_by_group.tsv (the original script reads
## "Supplementary Table_v6.xlsx", which was renamed). Same grammar and statistics,
## re-rendered narrower so panel a can be widened without downscaling this one.
suppressMessages({library(data.table); library(ggplot2)})
D <- "/Volumes/S840/04mypaper/lin/conflictresolution"
d <- fread(file.path(D,"fusio2trait/fusios_full_20260824/tables_fixed/recomb_by_group.tsv"))
LAB  <- c(opp="Opposite\nbetween sexes", conc="Concordant\nbetween sexes", bg="All tested\ngenes")
COLS <- setNames(c("#C0392B", "#3B6EA5", "#9A9A9A"), LAB)
d[, cls := factor(grp, levels=rev(LAB))]
p_oc <- wilcox.test(d[grp==LAB[["opp"]], rate], d[grp==LAB[["conc"]], rate])$p.value
p_ob <- wilcox.test(d[grp==LAB[["opp"]], rate], d[grp==LAB[["bg"]], rate])$p.value
cat(sprintf("opposite vs concordant P=%.2g | vs background P=%.2g\n", p_oc, p_ob))
pmP <- function(p){ e <- floor(log10(p)); m <- signif(p/10^e, 2)
  sprintf("italic(P)==%s%%*%%10^%d", m, e) }
med <- d[, .(m=median(rate)), by=cls]
set.seed(1)
p <- ggplot(d, aes(pmax(rate, 0.02), cls)) +
  geom_violin(aes(fill=cls), colour=NA, alpha=0.22, scale="width", width=0.58) +
  geom_jitter(data=d[cls!=LAB[["bg"]]], aes(colour=cls), height=0.10, size=0.55, alpha=0.30, stroke=0) +
  geom_boxplot(width=0.10, outlier.shape=NA, linewidth=0.34, fill="white", colour="grey30") +
  geom_point(data=med, aes(m, cls), shape=23, size=1.7, fill="white", colour="black", stroke=0.55) +
  annotate("segment", x=5.2, xend=5.2, y=2, yend=3, linewidth=0.3, colour="grey35") +
  annotate("segment", x=4.2, xend=5.2, y=2, yend=2, linewidth=0.3, colour="grey35") +
  annotate("segment", x=4.2, xend=5.2, y=3, yend=3, linewidth=0.3, colour="grey35") +
  annotate("text", x=6.6, y=2.5, label=pmP(p_oc), parse=TRUE, size=2.2, angle=90, colour="grey20") +
  annotate("segment", x=11.5, xend=11.5, y=1, yend=3, linewidth=0.3, colour="grey35") +
  annotate("segment", x=9.3, xend=11.5, y=1, yend=1, linewidth=0.3, colour="grey35") +
  annotate("segment", x=9.3, xend=11.5, y=3, yend=3, linewidth=0.3, colour="grey35") +
  annotate("text", x=14.8, y=2, label=pmP(p_ob), parse=TRUE, size=2.2, angle=90, colour="grey20") +
  scale_fill_manual(values=COLS, guide="none") + scale_colour_manual(values=COLS, guide="none") +
  scale_x_log10(name="local recombination rate (cM/Mb)",
                breaks=c(0.03,0.1,0.3,1,3,10), labels=c("0.03","0.1","0.3","1","3","10")) +
  labs(y=NULL) + coord_cartesian(xlim=c(0.02, 26), ylim=c(0.45, 3.9), clip="on") +
  theme_classic(base_size=7.4) +
  theme(axis.text.y=element_text(size=6.6, lineheight=0.95),
        axis.title.x=element_text(size=7.2), plot.margin=margin(4,6,3,3))
fix_panel <- function(p, w_pt, h_pt=150){ g <- ggplotGrob(p)
  g$widths[g$layout$l[g$layout$name=="panel"]]  <- grid::unit(w_pt,"pt")
  g$heights[g$layout$t[g$layout$name=="panel"]] <- grid::unit(h_pt,"pt"); g }
dev <- if (capabilities("cairo")) cairo_pdf else pdf
ggsave(file.path(D,"fusio2trait/figures_final/Fig5b_recomb.pdf"), fix_panel(p, 119),
       width=205/72, height=226/72, device=dev, bg="white")
cat("written Fig5b_recomb.pdf at 205x226 pt\n")
