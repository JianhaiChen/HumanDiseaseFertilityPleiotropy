#!/usr/bin/env Rscript
## Rebuild the directional-class figures at Benjamini-Hochberg FDR < 0.05
## (previously nominal P < 0.05). Recovered from the session transcripts and
## re-pointed at mode=="fdr", q==0.05 in antag50_classes.tsv.
##   Fig2D_antagclass.pdf      -> main Fig 4a
##   FigS_classprev.pdf        -> Supplementary Fig 10
##   FigS_finngen_rg4class.pdf -> Supplementary Fig 13
suppressMessages({library(data.table); library(ggplot2); library(cowplot)})
D  <- "/Volumes/S840/04mypaper/lin/conflictresolution"
CL <- file.path(D,"fusio2trait/sv_antisense_analysis_20260813/antag50_classes.tsv")
OUT<- file.path(D,"fusio2trait/figures_final")
COLS <- c(Male="#2C5985", Female="#C0392B")
lv <- c("d+f+","d-f-","d+f-","d-f+")
pmP <- function(p){p<-signif(p,2); e<-floor(log10(p)+1e-9)
  ifelse(p<1e-3, sprintf("%.1f%%*%%10^%d", p/10^e, e), sprintf("'%.2g'", p))}
dev <- if (capabilities("cairo")) cairo_pdf else pdf

cl <- fread(CL)[mode=="fdr" & q==0.05]
st <- fread(file.path(D,"fig5_data/section_table.tsv"))[, .(phenocode, prev)][!is.na(prev) & prev>0]
m  <- merge(cl, st, by="phenocode"); m[, lp := log10(prev)]

## ---- Fig 4a: pooled antagonistic share -------------------------------------
long <- rbind(m[n_m>=20, .(lp, axis="Male",   sh=100*(pp_m+mm_m)/n_m)],
              m[n_f>=20, .(lp, axis="Female", sh=100*(pp_f+mm_f)/n_f)])
sB <- long[, {ct<-cor.test(sh,lp); .(r=ct$estimate,p=ct$p.value,n=.N)}, by=axis]
print(sB)
sB[, lab := sprintf("'%s: r=%.2f, '*italic(P)==%s", substr(axis,1,1), r, pmP(p))]
pD <- ggplot(long, aes(sh, lp, colour=axis)) +
  geom_vline(xintercept=50, colour="grey85", linetype=2, linewidth=0.3) +
  geom_point(size=0.7, alpha=0.35, stroke=0) +
  geom_smooth(aes(linetype=axis), method="lm", formula=y~x, se=FALSE, linewidth=0.8) +
  geom_text(data=sB[axis=="Male"],   aes(x=-Inf,y=Inf,label=lab), hjust=-0.05, vjust=1.5, size=2.4, parse=TRUE, show.legend=FALSE) +
  geom_text(data=sB[axis=="Female"], aes(x=-Inf,y=Inf,label=lab), hjust=-0.05, vjust=3,   size=2.4, parse=TRUE, show.legend=FALSE) +
  scale_colour_manual(values=COLS, name=NULL, labels=c(Male="Male fertility", Female="Female fertility")) +
  scale_linetype_manual(values=c(Male="solid", Female="22"), guide="none") +
  scale_y_continuous(breaks=log10(c(0.01,0.1,1,10)), labels=c("0.01","0.1","1","10"), name="Population prevalence (%)") +
  scale_x_continuous(name="disease-fertility antagonistic pleiotropic genes (%)") +
  guides(colour=guide_legend(override.aes=list(size=2, alpha=1))) +
  theme_classic(base_size=10) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=9),
        legend.position="inside", legend.position.inside=c(0.80,0.14),
        legend.text=element_text(size=7.5), legend.key.size=unit(3.5,"mm"), legend.background=element_blank())
ggsave(file.path(OUT,"Fig2D_antagclass.pdf"), pD, width=89, height=80, units="mm", device=dev, bg="white")

## ---- Supplementary Fig 10: the four classes resolved ------------------------
cls <- rbind(
  m[n_m>=20, .(lp, axis="Male",   `d+f+`=100*pp_m/n_m, `d-f-`=100*mm_m/n_m,
               `d+f-`=100*dpfm_m/n_m, `d-f+`=100*dmfp_m/n_m)],
  m[n_f>=20, .(lp, axis="Female", `d+f+`=100*pp_f/n_f, `d-f-`=100*mm_f/n_f,
               `d+f-`=100*dpfm_f/n_f, `d-f+`=100*dmfp_f/n_f)])
lg <- melt(cls, id.vars=c("lp","axis"), variable.name="class", value.name="sh")
lg[, class := factor(class, levels=lv)]
sS <- lg[, {ct<-cor.test(sh,lp); .(r=ct$estimate,p=ct$p.value)}, by=.(class,axis)]
print(sS[order(class,axis)])
sS[, lab := sprintf("'%s: r=%+.2f, '*italic(P)==%s", substr(axis,1,1), r, pmP(p))]
sS[, vj := fifelse(axis=="Male", 1.4, 2.9)]
pS <- ggplot(lg, aes(sh, lp, colour=axis)) +
  geom_point(size=0.5, alpha=0.3, stroke=0) +
  geom_smooth(aes(linetype=axis), method="lm", formula=y~x, se=FALSE, linewidth=0.7) +
  geom_text(data=sS, aes(x=Inf,y=Inf,label=lab,vjust=vj), hjust=1.05, size=2.3, parse=TRUE, show.legend=FALSE) +
  facet_wrap(~class, nrow=2, scales="free_x") +
  scale_colour_manual(values=COLS, name=NULL, labels=c(Male="Male fertility", Female="Female fertility")) +
  scale_linetype_manual(values=c(Male="solid", Female="22"), guide="none") +
  scale_y_continuous(breaks=log10(c(0.01,0.1,1,10)), labels=c("0.01","0.1","1","10"), name="Population prevalence (%)") +
  scale_x_continuous(name="class share among shared genes (%)") +
  theme_bw(base_size=9) +
  theme(panel.grid=element_blank(), strip.background=element_blank(),
        strip.text=element_text(size=8.5, face="bold"),
        axis.text=element_text(size=7), axis.title=element_text(size=8),
        legend.position="top", legend.text=element_text(size=7.5),
        legend.key.size=unit(3,"mm"), legend.margin=margin(0,0,-4,0))
ggsave(file.path(OUT,"FigS_classprev.pdf"), pS, width=160, height=150, units="mm", device=dev, bg="white")

## ---- Supplementary Fig 13: class share vs coupling --------------------------
x  <- fread(CL)[mode=="fdr" & q==0.05]
st2<- fread(file.path(D,"fig5_data/section_table.tsv"))[, .(phenocode, rhoM=rho_father, rhoF=rho_mother, rgM=rg_father, rgF=rg_mother)]
x  <- merge(x, st2, by="phenocode")
mk <- function(sx){ n <- x[[paste0("n_",sx)]]; keep <- n>=20
  cc <- if (sx=="m") x$rhoM else x$rhoF; gg <- if (sx=="m") x$rgM else x$rgF
  rbindlist(lapply(1:4, function(k){
    cnt <- x[[c(paste0("pp_",sx),paste0("mm_",sx),paste0("dpfm_",sx),paste0("dmfp_",sx))[k]]]
    data.table(axis=ifelse(sx=="m","Male","Female"), class=lv[k],
               sh=100*cnt[keep]/n[keep], rho=cc[keep], rg=gg[keep])}))}
d <- rbind(mk("m"), mk("f")); d[, class := factor(class, levels=lv)]
lng <- rbind(d[,.(axis,class,sh,x=rg,meas="A_rg")][is.finite(x)],
             d[,.(axis,class,sh,x=rho,meas="B_rho")][is.finite(x)])
stats <- lng[, {ct<-cor.test(sh,x); .(r=ct$estimate,p=ct$p.value,n=.N)}, by=.(meas,class,axis)]
print(stats[order(meas,class,axis)], nrows=40)
stats[, lab := sprintf("'%s: r=%+.2f, '*italic(P)==%s", substr(axis,1,1), r, pmP(p))]
stats[, vj := fifelse(axis=="Male", 1.4, 2.9)]
labx <- c(A_rg=expression(italic(r)[g]*" (LDSC, genome-wide SNPs)"),
          B_rho=expression(rho[GE]*" (MR, expression effect)"))
for (mm in c("A_rg","B_rho")) {
  pp <- ggplot(lng[meas==mm], aes(x, sh, colour=axis)) +
    geom_point(size=0.35, alpha=0.18, stroke=0) +
    geom_smooth(aes(linetype=axis), method="lm", formula=y~x, se=FALSE, linewidth=0.7) +
    geom_text(data=stats[meas==mm], aes(x=Inf,y=Inf,label=lab,vjust=vj), hjust=1.05, size=2.2, parse=TRUE, show.legend=FALSE) +
    facet_wrap(~class, nrow=1, scales="free_x") +
    scale_colour_manual(values=COLS, name=NULL, labels=c(Male="Male fertility", Female="Female fertility")) +
    scale_linetype_manual(values=c(Male="solid", Female="22"), guide="none") +
    labs(x=labx[[mm]], y="class share among\nshared genes (%)") +
    theme_bw(base_size=9) +
    theme(panel.grid=element_blank(), strip.background=element_blank(),
          strip.text=element_text(size=8.5, face="bold"),
          axis.text=element_text(size=7), axis.title=element_text(size=8),
          legend.position="top", legend.text=element_text(size=7.5),
          legend.key.size=unit(3,"mm"), legend.margin=margin(0,0,-4,0))
  assign(paste0("p_",mm), pp)
}
fig <- plot_grid(p_A_rg, p_B_rho, ncol=1, labels=c("a","b"), label_size=11)
ggsave(file.path(OUT,"FigS_finngen_rg4class.pdf"), fig, width=180, height=150, units="mm", device=dev, bg="white")
cat("\nall three rebuilt at FDR 0.05\n")
