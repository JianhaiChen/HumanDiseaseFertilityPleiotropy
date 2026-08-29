#!/usr/bin/env Rscript
## Fig 4 panels B (H_g heterogeneity) and C (OEP ~ coupling concordance),
## recomputed from the manuscript's described logic and drawn in house style.
## Validated against the collaborator figure: H_g 302/1,828, 0.811/0.615
## (published 303/1,831, 0.810/0.614); OEP 1,542 pairs / 80 diseases,
## weighted partial r = -0.24 (published -0.23).
suppressMessages({library(data.table); library(ggplot2)})
MR <- "${MR}"
D  <- "${ROOT}"
S  <- "${SCRATCH}"
OUT <- file.path(D, "fusio2trait/figures_final")
strip <- function(x) sub("[.].*","",x)
pmP <- function(p){ p <- signif(p, 2); e <- floor(log10(p) + 1e-9)
  ifelse(p < 1e-3, sprintf("%.1f%%*%%10^%d", p/10^e, e), sprintf("'%.2g'", p)) }

## ================= panel B: H_g =================
b1g <- fread(file.path(MR,"gene_level_b1_disease_acat_summary.dropbox.csv"),
             select=c("gene","disease","p_bonf_global","top_beta"))
b1g[, gene := strip(gene)]
info <- fread(file.path(MR,"disease_gwas_source_info.csv"))[, .(disease, cat=Category)]
b2 <- fread(file.path(MR,"gene_level_b2_fitness_acat_summary.dropbox.csv"),
            select=c("gene","fitness","p_bonf_global","top_beta","top_se"))
b2[, gene := strip(gene)]
fert <- unique(b2[p_bonf_global<0.05, gene])
x <- merge(b1g[p_bonf_global<0.05 & is.finite(top_beta)], info, by="disease")
g <- x[, .(n=.N, p=mean(top_beta>0)), by=gene][n>=2]
win <- x[, .(nc=.N, pc=mean(top_beta>0)), by=.(gene,cat)]
g <- merge(g, win[, .(Hw = sum(nc*4*pc*(1-pc))/sum(nc)), by=gene], by="gene")
g[, Ht := 4*p*(1-p)]; g[, Hb := Ht - Hw]
g[, grp := fifelse(gene %in% fert, "Fertility-associated", "Other genes")]
sm <- rbindlist(lapply(c(Total="Ht", `Within category`="Hw", `Between category`="Hb"),
  function(v) g[, .(m=mean(get(v)), se=sd(get(v))/sqrt(.N)), by=grp]), idcol="comp")
pvB <- sapply(c("Ht","Hw","Hb"), function(v) summary(lm(get(v) ~ grp, g))$coefficients[2,4])
cat("panel B:", g[grp=="Fertility-associated",.N], "vs", g[grp=="Other genes",.N], "\n")
print(sm); cat("P:", signif(pvB,2), "\n")
sm[, comp := factor(comp, levels=c("Total","Within category","Between category"))]
ann <- data.table(comp=factor(c("Total","Within category","Between category"),
                              levels=levels(sm$comp)),
                  y = sm[, max(m+1.96*se), by=comp]$V1 + 0.055,
                  lab = sprintf("italic(P)==%s", pmP(pvB)))
pB2 <- ggplot(sm, aes(comp, m, fill=grp)) +
  geom_col(position=position_dodge(0.7), width=0.62, alpha=0.9) +
  geom_errorbar(aes(ymin=m-1.96*se, ymax=m+1.96*se), width=0.15,
                position=position_dodge(0.7), linewidth=0.4, colour="grey25") +
  geom_text(data=ann, aes(comp, y, label=lab), inherit.aes=FALSE, parse=TRUE, size=2.1) +
  scale_x_discrete(labels=c("Total","Within\ncategory","Between\ncategory")) +
  scale_fill_manual(values=c(`Fertility-associated`="#2C5985", `Other genes`="grey65"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08)),
                     name=expression("directional heterogeneity ("*italic(H)[g]*")")) +
  labs(x=NULL) +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(size=7.5), axis.text.y=element_text(size=8),
        axis.title=element_text(size=9), legend.position="top",
        legend.text=element_text(size=7), legend.key.size=unit(3,"mm"),
        legend.margin=margin(0,0,-4,0))
ggsave(file.path(OUT,"Fig4_panel_Hg.pdf"), pB2, width=89, height=80, units="mm",
       device=cairo_pdf, bg="white")

## ================= panel C: OEP ~ concordance =================
co <- fread(file.path(D,"fig5_data/gene_coords_v19.tsv"), header=FALSE,
            col.names=c("gene","chr","pos")); co[, gene := strip(gene)]
mhc  <- co[chr=="6" & pos>25e6 & pos<34e6, unique(gene)]
auto <- co[chr %in% as.character(1:22), unique(gene)]
b1 <- fread(file.path(MR,"gene_level_b1_disease_acat_summary.dropbox.csv"),
            select=c("gene","disease","p_bonf_across_genes","top_beta","top_se"))
b1[, gene := strip(gene)]
b1 <- b1[is.finite(top_beta) & is.finite(top_se) & top_se>0 & gene %in% auto & !gene %in% mhc]
b1[, z := top_beta/top_se]
norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
info2 <- fread(file.path(MR,"disease_gwas_source_info.csv"))[, .(disease, key=norm(Disease))]
fix <- data.table(kk = norm(c("Blood Clot Lung","Lung Cancer","Renal Cancer","Esophageal Cancer",
        "Colorectal Cancer","Prostate Cancer","Ovarian Cancer","Polycystic Ovary Syndrome",
        "Breast Cancer","Endometrial Cancer")),
  disease = c("Blood_Clot_Lung","lungcarcinoma","renalcancerwu","esophagealcancerwu",
        "colorectalcancerwu","GCST90274714Prostatecancer","ovariancancerwu",
        "Polycystic_ovary_syndrome","EFO_0000305Breastcancer.h.tsv.gz1","endometrialcancerwu"))
setnames(fix,"kk","key"); info2 <- rbind(info2[!key %in% fix$key], fix)
s3 <- fread(file.path(S,"tableS3.tsv"))[, .(key=norm(Disease), rhom=rhoGE_male, rhof=rhoGE_female)]
s3 <- merge(s3, info2, by="key")
sig <- b1[p_bonf_across_genes < 0.05]
ds <- split(sig[, .(gene, s = sign(z))], sig$disease); dz <- names(ds)
pairs <- list()
for (i in seq_along(dz)) for (j in seq_len(i-1)) {
  m <- merge(ds[[dz[i]]], ds[[dz[j]]], by="gene")
  if (nrow(m) >= 30) pairs[[length(pairs)+1]] <-
    data.table(d1=dz[i], d2=dz[j], n=nrow(m), oep=mean(m$s.x != m$s.y))
}
P <- rbindlist(pairs)
P <- merge(P, s3[,.(d1=disease, r1m=rhom, r1f=rhof)], by="d1")
P <- merge(P, s3[,.(d2=disease, r2m=rhom, r2f=rhof)], by="d2")
P[, q := (r1m*r2m + r1f*r2f)/2]
Z <- dcast(b1[abs(z)<9], gene ~ disease, value.var="z")
CC <- cor(as.matrix(Z[,-1]), use="pairwise.complete.obs")
P[, Sij := CC[cbind(d1,d2)]]
P <- P[is.finite(q) & is.finite(Sij)]
wres <- function(y, x, w) resid(lm(y ~ x, weights=w))
P[, `:=`(ry = wres(oep, scale(Sij)[,1], n), rx = wres(scale(q)[,1], scale(Sij)[,1], n))]
wpr <- cov.wt(P[,cbind(ry,rx)], wt=P$n, cor=TRUE)$cor[1,2]
## permutation P (disease-label permutation of rho)
set.seed(1); B <- 100000
rhos <- s3[disease %in% unique(c(P$d1,P$d2))]
obs <- abs(wpr); hits <- 0
d1i <- match(P$d1, rhos$disease); d2i <- match(P$d2, rhos$disease)
Ss <- scale(P$Sij)[,1]
for (b in 1:B) {
  pi <- sample(nrow(rhos))
  qb <- (rhos$rhom[pi][d1i]*rhos$rhom[pi][d2i] + rhos$rhof[pi][d1i]*rhos$rhof[pi][d2i])/2
  rxb <- wres(scale(qb)[,1], Ss, P$n)
  if (abs(cov.wt(cbind(P$ry, rxb), wt=P$n, cor=TRUE)$cor[1,2]) >= obs) hits <- hits+1
}
permP <- (1+hits)/(B+1)
cat(sprintf("panel C: %d pairs / %d diseases; weighted partial r=%+.3f; perm P=%.2g\n",
    nrow(P), uniqueN(c(P$d1,P$d2)), wpr, permP))

## decile-binned residual plot
P[, dec := cut(rx, quantile(rx, 0:10/10), include.lowest=TRUE, labels=FALSE)]
bins <- P[, .(x = weighted.mean(rx, n), y = weighted.mean(ry, n),
              se = sqrt(sum(n*(ry-weighted.mean(ry,n))^2)/sum(n))/sqrt(.N)), by=dec]
fit <- lm(ry ~ rx, data=P, weights=P$n)
labC <- sprintf("'weighted partial r=%.2f, '*italic(P)==%s", wpr, pmP(permP))
pC2 <- ggplot(bins, aes(x, y)) +
  geom_hline(yintercept=0, colour="grey85", linewidth=0.3) +
  geom_abline(intercept=coef(fit)[1], slope=coef(fit)[2], colour="#C0392B", linewidth=0.7) +
  geom_errorbar(aes(ymin=y-se, ymax=y+se), width=0, colour="#2C5985", linewidth=0.4) +
  geom_point(size=1.6, colour="#2C5985") +
  annotate("text", x=-Inf, y=-Inf, hjust=-0.05, vjust=-0.8, size=2.3, parse=TRUE, label=labC) +
  annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=1.5, size=2.2, colour="grey40",
           label=sprintf("%d pairs / %d diseases", nrow(P), uniqueN(c(P$d1,P$d2)))) +
  labs(x="residual coupling concordance", y="residual opposite-effect proportion") +
  theme_classic(base_size=10) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8.5))
ggsave(file.path(OUT,"Fig4_panel_OEP.pdf"), pC2, width=89, height=80, units="mm",
       device=cairo_pdf, bg="white")
cat("panels written\n")
