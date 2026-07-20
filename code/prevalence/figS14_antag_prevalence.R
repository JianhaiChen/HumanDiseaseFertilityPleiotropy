## SUPP figure S14: gene-level %antagonistic pleiotropic genes vs disease prevalence (FinnGen, n=2149).
## Antagonistic = concordant-sign (d+f+/d-f-) genes significant (nominal P<0.05) for both disease & fertility;
## %antag per endpoint correlated with log10 FinRegistry prevalence. Panels C=Male axis, D=Female axis.
#!/usr/bin/env Rscript
# Fig5 C/D: gene-level %antagonistic pleiotropy vs FinnGen disease prevalence.
# antag = concordant sign (d+f+ / d-f-); syn = discordant. pct_antag = 100*antag/(antag+syn).
# Hypothesis: antagonistic (trade-off) pleiotropy tracks HIGHER prevalence; synergistic lower.
suppressMessages(library(data.table))
D  <- "/Volumes/S840/04mypaper/lin/conflictresolution/rg_prevalence_analysis_2026-07-16"
an <- fread(file.path(D,"tables","antag_syn_ALL.tsv"))            # fg axis thr n_both n_antag n_syn pct_antag
pv <- fread(file.path(D,"tables","finregistry_prevalence_R12.tsv"))
pv <- pv[finregistry_available==TRUE,
         .(fg=endpoint_id, prev=as.numeric(whole_population_prevalence_all))]
pv <- pv[is.finite(prev) & prev>0][, lp := log10(prev)]

FLOOR <- 30          # min genes-in-both for a stable proportion (revisit from table below)
m <- merge(an, pv, by="fg")

# ---- scan: r(pct_antag, log10 prev) over axis x threshold x floor ----
cat("axis  thr     floor  n_dz   r      P\n")
for(ax in c("father","mother")) for(th in sort(unique(an$thr))) for(fl in c(10,20,30,50)){
  s <- m[axis==ax & thr==th & n_both>=fl & is.finite(pct_antag)]
  if(nrow(s)<20) next
  ct <- cor.test(s$pct_antag, s$lp)
  cat(sprintf("%-6s %-7g %-5d %5d  %+.3f  %.2e\n", ax, th, fl, nrow(s), ct$estimate, ct$p.value))
}

# ---- pick main: thr=0.05, floor=30 (largest per-disease gene set) ----
TH <- 0.05
mk <- function(ax){
  s <- m[axis==ax & thr==TH & n_both>=FLOOR & is.finite(pct_antag)]
  s[, wt := n_both]
  s
}
fa <- mk("father"); mo <- mk("mother")
saveRDS(list(fa=fa,mo=mo), file.path(D,"antag_prev_panels.rds"))

panel <- function(s, title, col){
  ct <- cor.test(s$pct_antag, s$lp)
  fit <- lm(pct_antag ~ lp, data=s)
  # prevalence tertiles median (house style)
  s[, tb := cut(lp, quantile(lp, c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("low","mid","high"))]
  med <- s[, .(x=median(lp), y=median(pct_antag)), by=tb][order(x)]
  xr <- range(s$lp); yr <- c(min(30,quantile(s$pct_antag,.01)), max(70,quantile(s$pct_antag,.99)))
  plot(s$lp, s$pct_antag, pch=16, col=adjustcolor(col,.28),
       cex=0.5+1.2*sqrt(s$n_both/max(s$n_both)),
       xlab=expression(log[10]~"disease prevalence (%)"), ylab="% antagonistic genes",
       main=title, ylim=yr, las=1, bty="l")
  abline(h=50, lty=3, col="grey60")
  xs <- seq(xr[1],xr[2],length=100); pr <- predict(fit, data.frame(lp=xs), interval="confidence")
  polygon(c(xs,rev(xs)), c(pr[,"lwr"],rev(pr[,"upr"])), col=adjustcolor(col,.15), border=NA)
  lines(xs, pr[,"fit"], col=col, lwd=2.5)
  points(med$x, med$y, pch=23, bg="white", col=col, cex=1.6, lwd=2)
  legend("topleft", bty="n", cex=0.95,
         legend=sprintf("r = %+.2f\nP = %.1e\nn = %d dz", ct$estimate, ct$p.value, nrow(s)))
}

pdf(file.path(D,"Fig5CD_antag_prevalence.pdf"), width=10, height=4.6)
par(mfrow=c(1,2), mar=c(4.5,4.5,3,1), cex.axis=0.95, cex.lab=1.05)
panel(fa, "C  Male-fitness axis", "#2C6FB0")   # C (stronger)
panel(mo, "D  Female-fitness axis", "#C0392B")   # D
dev.off()

# also PNG for quick open
png(file.path(D,"Fig5CD_antag_prevalence.png"), width=1600, height=740, res=150)
par(mfrow=c(1,2), mar=c(4.5,4.5,3,1), cex.axis=0.95, cex.lab=1.05)
panel(fa, "C  Male-fitness axis", "#2C6FB0")
panel(mo, "D  Female-fitness axis", "#C0392B")
dev.off()
cat("\nWROTE Fig5CD_antag_prevalence.pdf/png ; n father=",nrow(fa)," mother=",nrow(mo),"\n")
