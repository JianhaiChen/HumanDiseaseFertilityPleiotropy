#!/usr/bin/env Rscript
## Fig 5b: local recombination rate around the cross-sex conflict genes.
## Rate = cM per Mb from the HapMap Phase II map (GRCh37, PLINK distribution),
## measured in a +/-500 kb window centred on each gene.
suppressMessages({library(data.table); library(ggplot2); library(openxlsx)})
D   <- "${ROOT}"
MR  <- "${MR}/recomb"
OUT <- file.path(D, "fusio2trait/figures_final")

co <- fread(file.path(D,"fig5_data/gene_coords_v19.tsv"), header=FALSE,
            col.names=c("gene","chr","pos"))[, gene := sub("[.].*","",gene)]
co[, chr := suppressWarnings(as.integer(chr))]
rec <- rbindlist(lapply(1:22, function(c) {
  r <- fread(file.path(MR, sprintf("plink.chr%d.GRCh37.map", c)), header=FALSE)
  data.table(chr=c, bp=r$V4, cm=r$V3) }))
rate_of <- function(dt) { dt[, rate := NA_real_]
  for (c in unique(dt$chr)) { m <- rec[chr==c][order(bp)]; i <- dt$chr==c
    dt$rate[i] <- approx(m$bp, m$cm, pmin(dt$pos[i]+5e5, max(m$bp)), rule=2)$y -
                  approx(m$bp, m$cm, pmax(dt$pos[i]-5e5, min(m$bp)), rule=2)$y }
  dt[is.finite(rate)] }

S9 <- setDT(read.xlsx(file.path(D,"fusio2trait/Supplementary Table_v6.xlsx"),
                      sheet="Supplementary Table 9", startRow=2))
S9[, `:=`(gene=sub("[.].*","",gene), opp=sign_class_metaEffect=="opposite")]
fe <- rbindlist(lapply(list.files(file.path(D,"fusio2trait/fusios_full_20260824/acat50_fert_full"),
                                  "tsv.gz$", full.names=TRUE),
                       function(f) fread(cmd=paste("gzcat", f), select="gene")))
tested <- unique(sub("[.].*","", fe$gene))

g  <- rate_of(merge(S9[, .(gene, opp)], co, by="gene")[!is.na(chr)])
set.seed(1)
bgg <- rate_of(co[gene %in% tested & !gene %in% g$gene & !is.na(chr)][sample(.N, 8000)])
d <- rbind(bgg[, .(rate, grp="All tested\ngenes")],
           g[opp==FALSE, .(rate, grp="Concordant\nbetween sexes")],
           g[opp==TRUE,  .(rate, grp="Opposite\nbetween sexes")])
d[, grp := factor(grp, levels=c("All tested\ngenes","Concordant\nbetween sexes","Opposite\nbetween sexes"))]
fwrite(d, file.path(D,"fusio2trait/fusios_full_20260824/tables_fixed/recomb_by_group.tsv"), sep="\t")

p1 <- wilcox.test(g[opp==TRUE, rate], g[opp==FALSE, rate])$p.value
p2 <- wilcox.test(g[opp==TRUE, rate], bgg$rate)$p.value
lab <- function(p) { e <- floor(log10(p))
  paste0(sprintf("%.0f", p/10^e), "×10", chartr("-0123456789","⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(e))) }
cat(sprintf("n: bg %d | concordant %d | opposite %d\n", nrow(bgg), g[opp==FALSE,.N], g[opp==TRUE,.N]))
print(d[, .(median=round(median(rate),3), pct_below_0.5=round(100*mean(rate<0.5))), by=grp])
cat(sprintf("opposite vs concordant P=%.2g | vs background P=%.2g\n", p1, p2))

## same visual grammar as panel d: horizontal violins, box, jittered points,
## white median diamond, one P label per row, theme_classic
LAB  <- c(opp="Opposite\nbetween sexes", conc="Concordant\nbetween sexes", bg="All tested\ngenes")
COLS <- setNames(c("#C0392B", "#3B6EA5", "#9A9A9A"), LAB)
d[, cls := factor(grp, levels=rev(LAB))]
## panel b compares groups with each other, so the statistics are drawn as
## pairwise brackets rather than as one label per row (panel d tests each row against zero)
p_oc <- wilcox.test(g[opp==TRUE, rate], g[opp==FALSE, rate])$p.value
p_ob <- wilcox.test(g[opp==TRUE, rate], bgg$rate)$p.value
pmP <- function(p) { e <- floor(log10(p)); m <- signif(p/10^e, 2)
  sprintf("italic(P)==%s%%*%%10^%d", m, e) }
med <- d[, .(m=median(rate)), by=cls]
set.seed(1)
p <- ggplot(d, aes(pmax(rate, 0.02), cls)) +
  geom_violin(aes(fill=cls), colour=NA, alpha=0.22, scale="width", width=0.58) +
  geom_jitter(data=d[cls!=LAB[["bg"]]], aes(colour=cls), height=0.10, size=0.55,
              alpha=0.30, stroke=0) +
  geom_boxplot(width=0.10, outlier.shape=NA, linewidth=0.34, fill="white", colour="grey30") +
  geom_point(data=med, aes(m, cls), shape=23, size=1.7, fill="white", colour="black", stroke=0.55) +
  ## brackets: opposite vs concordant, opposite vs all tested genes
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
  labs(y=NULL) +
  coord_cartesian(xlim=c(0.02, 26), ylim=c(0.45, 3.9), clip="on") +
  theme_classic(base_size=7.4) +
  theme(axis.text.y=element_text(size=6.6, lineheight=0.95),
        axis.title.x=element_text(size=7.2), plot.margin=margin(4,6,3,3))

## force both horizontal-violin panels to the same plotting width so their x axes
## are the same length regardless of how wide the y-axis labels are
fix_panel <- function(p, w_pt = 166, h_pt = 150) {
  g <- ggplotGrob(p)
  g$widths[ g$layout$l[g$layout$name == "panel"] ]  <- grid::unit(w_pt, "pt")
  g$heights[ g$layout$t[g$layout$name == "panel"] ] <- grid::unit(h_pt, "pt")
  g
}

ggsave(file.path(OUT,"Fig5b_recomb.pdf"), fix_panel(p), width=252/72, height=226/72, device=cairo_pdf, bg="white")
cat("written Fig5b_recomb.pdf (panel-d style)\n")
