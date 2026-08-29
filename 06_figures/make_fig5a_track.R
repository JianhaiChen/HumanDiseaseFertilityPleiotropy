#!/usr/bin/env Rscript
## Fig 5a: two horizontal tracks, male fertility on top and female below.
## x = cross-tissue inverse-variance meta effect on fertility; a gene whose
## effect changes sign between the sexes crosses the centre line.
## Gene selection follows the volcano panel (assemble_fig4_fig5.R): non-MHC
## genes with a finite IVW meta effect in both sexes, BH q<0.01 in each sex.
suppressMessages({library(data.table); library(ggplot2)})
BASE <- "${ROOT}"
OUT  <- file.path(BASE,"fusio2trait/figures_final")
fix_panel <- function(p, w=166, h=150) { g <- ggplotGrob(p)
  g$widths[g$layout$l[g$layout$name=="panel"]] <- grid::unit(w,"pt")
  g$heights[g$layout$t[g$layout$name=="panel"]] <- grid::unit(h,"pt"); g }

d  <- fread(file.path(BASE,"fig5_data/b2_meta.csv"))
co <- fread(file.path(BASE,"fig5_data/gene_coords_v19.tsv"), header=FALSE,
            col.names=c("gene","chr","start"))
d <- merge(d, co, by="gene", all.x=TRUE)
d[, mhc := !is.na(chr) & chr==6 & start>=25e6 & start<=34e6]
x <- d[mhc==FALSE & is.finite(p_acat)]
x[, dir := fifelse(is.finite(meta_beta_ivw) & meta_beta_ivw!=0, sign(meta_beta_ivw), sign(top_beta))]
w <- dcast(x, gene ~ fitness, value.var=c("p_acat","dir","meta_beta_ivw"))
w <- w[is.finite(meta_beta_ivw_male) & is.finite(meta_beta_ivw_female)]
w[, `:=`(q_male=p.adjust(p_acat_male,"BH"), q_female=p.adjust(p_acat_female,"BH"))]
sel <- w[q_male<0.01 & q_female<0.01]
sel[, opp := dir_male != dir_female]
cat(sprintf("both-significant %d | opposite %d (%.1f%%) | same %d\n",
            nrow(sel), sum(sel$opp), 100*mean(sel$opp), sum(!sel$opp)))

seg <- sel[order(opp), .(opp, xm=meta_beta_ivw_male, xf=meta_beta_ivw_female)]
pts <- rbind(seg[, .(opp, x=xm, y=1)], seg[, .(opp, x=xf, y=0)])
COL <- c(`TRUE`="#C0392B", `FALSE`="grey72")

RX <- range(c(seg$xm, seg$xf)); PAD <- 0.03*diff(RX)
p <- ggplot() +
  annotate("segment", x=RX[1]-PAD, xend=RX[2]+PAD, y=c(0,1), yend=c(0,1),
           colour="grey60", linewidth=0.4) +
  geom_segment(data=seg, aes(x=xm, xend=xf, y=1, yend=0, colour=opp),
               linewidth=0.25, alpha=0.75) +
  geom_point(data=pts, aes(x, y, colour=opp), size=0.7, stroke=0) +
  annotate("segment", x=0, xend=0, y=-0.04, yend=1.04, colour="grey25", linewidth=0.35) +
  annotate("text", x=-Inf, y=1.24, hjust=-0.04, size=2.4, colour="#C0392B",
           label=sprintf("opposite sign: %d of %d (%.1f%%)",
                         sum(sel$opp), nrow(sel), 100*mean(sel$opp))) +
  annotate("text", x=-Inf, y=1.13, hjust=-0.04, size=2.4, colour="grey55",
           label=sprintf("same sign: %d (%.1f%%)", sum(!sel$opp), 100*mean(!sel$opp))) +
  scale_colour_manual(values=COL, guide="none") +
  scale_y_continuous(breaks=c(0,1), labels=c("Female\nfertility","Male\nfertility"),
                     limits=c(-0.10, 1.34), expand=c(0,0)) +
  scale_x_continuous(name=expression("cross-tissue meta effect on fertility ("*hat(beta)*")")) +
  labs(y=NULL) + theme_classic(base_size=7.4) +
  theme(axis.text.y=element_text(size=7, colour="black", lineheight=0.95),
        axis.text.x=element_text(size=6.8, colour="black"),
        axis.title.x=element_text(size=7.2), axis.line.y=element_blank(),
        axis.ticks.y=element_blank(), plot.margin=margin(6,6,3,3))
ggsave(file.path(OUT,"Fig5a_track.pdf"), fix_panel(p), width=252/72, height=226/72,
       device=cairo_pdf, bg="white")
cat("written Fig5a_track.pdf\n")
