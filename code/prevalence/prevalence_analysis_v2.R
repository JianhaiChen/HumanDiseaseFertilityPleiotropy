#!/usr/bin/env Rscript
# Disease-level rg x prevalence analysis (v2), pre-specified.
# Principle: mapping quality = VALIDITY (handled by exclusion/stratification/sensitivity);
#            GBD uncertainty = PRECISION (handled by weighting / Monte Carlo).
# No arbitrary mapping-quality weights in any inferential model.
suppressWarnings(suppressMessages({library(sandwich); library(lmtest); library(MASS)}))
set.seed(20260716)

BASE <- "/Volumes/S840/04mypaper/lin/conflictresolution/rg_prevalence_analysis_2026-07-16"
TAB  <- file.path(BASE, "tables"); FIG <- file.path(BASE, "figures")
rd <- function(f) read.delim(file.path(TAB, f), stringsAsFactors=FALSE, check.names=FALSE)

## ---- 1. Load ----
rg   <- rd("sexspec_rg_results.tsv")                              # endpoint,disease_sex,fitness,rg,se,z,p,h2..
prov <- rd("SuppTable_GWAS_GBD_prevalence_provenance.tsv")        # phenotype x sex: cause, mapping_type, prev mean/lo/hi
mast <- rd("MASTER_prevalence_table_for_GPT.tsv")                 # per GBD_cause: category, duration_flag
Nfg  <- tryCatch(read.delim(file.path(TAB,"fg_rg_jobs.tsv"), header=FALSE, stringsAsFactors=FALSE), error=function(e) NULL)
if(!is.null(Nfg)) names(Nfg) <- c("endpoint","N_eff")

## ---- 2. sex-matched rg per (endpoint, sex): male-disease x father ; female-disease x mother ----
rg$rg <- suppressWarnings(as.numeric(rg$rg)); rg$se <- suppressWarnings(as.numeric(rg$se))
mkey <- function(df, dsex, fit) {
  s <- df[df$disease_sex==dsex & df$fitness==fit, c("endpoint","rg","se")]
  names(s) <- c("phenotype", paste0("rg_",dsex), paste0("rgse_",dsex)); s
}
rgM <- mkey(rg,"male","father"); rgF <- mkey(rg,"female","mother")
rgw <- merge(rgM, rgF, by="phenotype", all=TRUE)

## ---- 3. prevalence wide per phenotype (Male/Female) ----
prov$prevalence_mean  <- suppressWarnings(as.numeric(prov$prevalence_mean))
prov$prevalence_lower <- suppressWarnings(as.numeric(prov$prevalence_lower))
prov$prevalence_upper <- suppressWarnings(as.numeric(prov$prevalence_upper))
pv <- prov[, c("GWAS_phenotype","GBD_cause_id","GBD_cause_name","mapping_type","GWAS_ICD10",
               "sex","prevalence_mean","prevalence_lower","prevalence_upper","exclusion_reason")]
names(pv)[1] <- "phenotype"
pw <- reshape(pv, idvar=c("phenotype","GBD_cause_id","GBD_cause_name","mapping_type","GWAS_ICD10","exclusion_reason"),
              timevar="sex", direction="wide")
# columns like prevalence_mean.Male etc.

## ---- 4. join rg + prevalence at phenotype level ----
d <- merge(pw, rgw, by="phenotype", all.x=TRUE)
d$source_type <- ifelse(grepl("^2000[12]_", d$phenotype), "self_report", "registry_icd")
d <- merge(d, if(!is.null(Nfg)) Nfg else data.frame(phenotype=character(),N_eff=numeric()),
           by.x="phenotype", by.y="endpoint", all.x=TRUE)

## keep only rows mapped to a GBD cause with a usable mapping_type
usable_map <- c("exact","narrower","broader")           # composite / none excluded from all_usable
d <- d[!is.na(d$GBD_cause_id) & d$GBD_cause_id!="" & d$mapping_type %in% usable_map, ]

## ---- 5. select ONE primary GWAS per analysis_disease_id (GBD_cause_id), outcome-independent ----
d$analysis_disease_id <- d$GBD_cause_id
prec  <- ifelse(d$mapping_type %in% c("exact","narrower"), 1L, 0L)   # rung2: precise vs broad
src   <- ifelse(d$source_type=="registry_icd", 1L, 0L)              # rung1: registry/ICD > self-report
gran  <- ifelse(d$mapping_type=="exact", 1L, 0L)                     # rung4: finer GBD granularity
Nv    <- ifelse(is.na(d$N_eff), -Inf, d$N_eff)                       # rung3: larger effective N
# sort key: src desc, prec desc, N desc, gran desc, name asc  (all pre-specified, no outcome terms)
ord <- order(d$analysis_disease_id, -src, -prec, -Nv, -gran, d$phenotype)
d <- d[ord, ]
d$is_primary <- !duplicated(d$analysis_disease_id)
prim <- d[d$is_primary, ]

## ---- 6. build long (disease x sex) analysis frame from primary GWAS ----
mk_long <- function(df) {
  do.call(rbind, lapply(c("Male","Female"), function(SX) {
    rgcol <- if(SX=="Male") "rg_male" else "rg_female"
    secol <- if(SX=="Male") "rgse_male" else "rgse_female"
    data.frame(
      analysis_disease_id = df$analysis_disease_id,
      GBD_cause_name = df$GBD_cause_name,
      phenotype = df$phenotype, source_type = df$source_type,
      mapping_type = df$mapping_type, sex = SX,
      rg = df[[rgcol]], rg_se = df[[secol]],
      prev = df[[paste0("prevalence_mean.",SX)]],
      prev_lo = df[[paste0("prevalence_lower.",SX)]],
      prev_hi = df[[paste0("prevalence_upper.",SX)]],
      stringsAsFactors = FALSE)
  }))
}
L <- mk_long(prim)

## disease metadata: category + duration/fatality from MASTER (validity strata, pre-specified)
m2 <- mast[, c("GBD_cause","category","duration_flag")]; names(m2)[1] <- "GBD_cause_name"
L <- merge(L, m2, by="GBD_cause_name", all.x=TRUE)
L$category[is.na(L$category)|L$category==""] <- "Other"
L$duration_flag[is.na(L$duration_flag)] <- ""
L$acute        <- L$duration_flag=="acute"
L$high_fatality <- L$duration_flag=="high-fatality"
L$is_infection <- L$category=="Infection"
L$is_aggregate <- L$category=="Other/Aggregate"
# coarse category for the minimally-adjusted model (n~15 cannot support 13 levels)
grp <- c(Cancer="Cancer", Cardiovascular="Cardiometabolic", "Metabolic/Renal"="Cardiometabolic",
         Neurological="Neuropsychiatric", Psychiatric="Neuropsychiatric")
L$cat_broad <- ifelse(L$category %in% names(grp), grp[L$category], "Other")

## outcome + prevalence-precision on log10 scale
L$log10prev <- log10(L$prev)
L$se_log10  <- (log10(L$prev_hi) - log10(L$prev_lo)) / (2*1.96)
L <- L[is.finite(L$rg) & is.finite(L$log10prev), ]

## ---- 7. datasets (pre-specified inclusion/exclusion) ----
# global validity screen: drop infections + aggregates from every dataset (not diseases of interest)
core <- L[!L$is_infection & !L$is_aggregate & is.finite(L$se_log10) & L$se_log10>0, ]
core$strict_exact       <- core$mapping_type=="exact"
core$exact_plus_close   <- core$mapping_type %in% c("exact","narrower")   # close := interpretable narrower
core$all_usable         <- core$mapping_type %in% c("exact","narrower","broader")
core$chronic_only       <- core$strict_exact & !core$acute
core$exclude_high_fatality <- core$strict_exact & !core$high_fatality
write.csv(core, file.path(TAB,"prevalence_analysis_datasets.csv"), row.names=FALSE)

## ---- helpers: HC3 regression summary ----
fit_hc3 <- function(dat, adj=FALSE, label, dataset, sex) {
  if(nrow(dat) < 5 || length(unique(dat$rg))<3) return(NULL)
  # minimally adjusted: coarse category + high_fatality; drop terms that are constant/singleton to avoid rank deficiency
  if(adj) {
    dat$cat_broad <- factor(dat$cat_broad)
    terms <- "rg"
    if(nlevels(droplevels(dat$cat_broad))>1 && min(table(dat$cat_broad))>=2) terms <- c(terms,"cat_broad")
    if(length(unique(dat$high_fatality))>1) terms <- c(terms,"high_fatality")
    f <- as.formula(paste("log10prev ~", paste(terms, collapse=" + ")))
  } else f <- log10prev ~ rg
  m <- lm(f, data=dat)
  V <- vcovHC(m, type="HC3"); ct <- coeftest(m, vcov.=V)
  b <- coef(m)["rg"]; se <- sqrt(V["rg","rg"]); df <- df.residual(m)
  tcrit <- qt(0.975, df)
  pear <- suppressWarnings(cor.test(dat$rg, dat$log10prev, method="pearson"))
  spear<- suppressWarnings(cor.test(dat$rg, dat$log10prev, method="spearman", exact=FALSE))
  data.frame(model=label, dataset=dataset, sex=sex, n=nrow(dat),
             beta=unname(b), se_hc3=unname(se),
             ci_lo=unname(b-tcrit*se), ci_hi=unname(b+tcrit*se),
             p=unname(ct["rg","Pr(>|t|)"]),
             pearson_r=unname(pear$estimate), pearson_p=pear$p.value,
             spearman_r=unname(spear$estimate), spearman_p=spear$p.value,
             stringsAsFactors=FALSE)
}
fit_rlm <- function(dat, label, dataset, sex) {
  if(nrow(dat) < 5 || length(unique(dat$rg))<3) return(NULL)
  m <- tryCatch(rlm(log10prev ~ rg, data=dat, maxit=100), error=function(e) NULL); if(is.null(m)) return(NULL)
  s <- summary(m)$coefficients; b <- s["rg","Value"]; se <- s["rg","Std. Error"]
  p <- 2*pnorm(-abs(b/se))
  data.frame(model=label, dataset=dataset, sex=sex, n=nrow(dat),
             beta=b, se_hc3=se, ci_lo=b-1.96*se, ci_hi=b+1.96*se, p=p,
             pearson_r=NA, pearson_p=NA, spearman_r=NA, spearman_p=NA, stringsAsFactors=FALSE)
}
# uncertainty-weighted (winsorized inverse-variance + tau2), sensitivity
fit_wls <- function(dat, label, dataset, sex) {
  if(nrow(dat) < 5 || length(unique(dat$rg))<3) return(NULL)
  tau2 <- median(dat$se_log10^2, na.rm=TRUE)
  w <- 1/(dat$se_log10^2 + tau2)
  qs <- quantile(w, c(.05,.95), na.rm=TRUE); w <- pmin(pmax(w, qs[1]), qs[2])  # winsorize
  m <- lm(log10prev ~ rg, data=dat, weights=w)
  V <- vcovHC(m, type="HC3"); b <- coef(m)["rg"]; se <- sqrt(V["rg","rg"]); df <- df.residual(m)
  tcrit <- qt(0.975, df)
  data.frame(model=label, dataset=dataset, sex=sex, n=nrow(dat),
             beta=unname(b), se_hc3=unname(se), ci_lo=unname(b-tcrit*se), ci_hi=unname(b+tcrit*se),
             p=unname(coeftest(m,vcov.=V)["rg","Pr(>|t|)"]),
             pearson_r=NA,pearson_p=NA,spearman_r=NA,spearman_p=NA, stringsAsFactors=FALSE)
}
sub <- function(flag, sx) core[core[[flag]] & core$sex==sx, ]

## ---- 8. PRIMARY: strict_exact, per sex, unadjusted + minimally adjusted ----
primary <- do.call(rbind, list(
  fit_hc3(sub("strict_exact","Male"),   FALSE, "unadjusted_HC3", "strict_exact","Male"),
  fit_hc3(sub("strict_exact","Female"), FALSE, "unadjusted_HC3", "strict_exact","Female"),
  fit_hc3(sub("strict_exact","Male"),   TRUE,  "minimally_adjusted_HC3", "strict_exact","Male"),
  fit_hc3(sub("strict_exact","Female"), TRUE,  "minimally_adjusted_HC3", "strict_exact","Female")
))
write.csv(primary, file.path(TAB,"prevalence_regression_primary.csv"), row.names=FALSE)

## ---- 9. SENSITIVITY: all datasets x {unadjusted, robust rlm, weighted} x sex + leave-one-mapping-out + pooled ----
DS <- c("strict_exact","exact_plus_close","all_usable","chronic_only","exclude_high_fatality")
sens <- list()
for(ds in DS) for(sx in c("Male","Female")) {
  dat <- sub(ds, sx)
  sens <- c(sens, list(fit_hc3(dat, FALSE, "unadjusted_HC3", ds, sx),
                       fit_rlm(dat, "robust_rlm", ds, sx),
                       fit_wls(dat, "uncertainty_weighted", ds, sx)))
}
# leave-one-mapping-class-out on all_usable
for(cls in c("exact","narrower","broader")) for(sx in c("Male","Female")) {
  dat <- core[core$all_usable & core$mapping_type!=cls & core$sex==sx, ]
  sens <- c(sens, list(fit_hc3(dat, FALSE, paste0("LOMO_drop_",cls), "all_usable", sx)))
}
# pooled both sexes, cluster-robust by disease (secondary)
poolfit <- function(dat, ds) {
  if(nrow(dat)<6) return(NULL)
  m <- lm(log10prev ~ rg + sex, data=dat)
  V <- vcovCL(m, cluster=~analysis_disease_id, type="HC1")
  b <- coef(m)["rg"]; se <- sqrt(V["rg","rg"]); df <- length(unique(dat$analysis_disease_id))-1
  tcrit <- qt(0.975, max(df,1))
  data.frame(model="pooled_clusterrobust", dataset=ds, sex="Both", n=nrow(dat),
             beta=unname(b), se_hc3=unname(se), ci_lo=unname(b-tcrit*se), ci_hi=unname(b+tcrit*se),
             p=unname(coeftest(m,vcov.=V)["rg","Pr(>|t|)"]),
             pearson_r=NA,pearson_p=NA,spearman_r=NA,spearman_p=NA, stringsAsFactors=FALSE)
}
sens <- c(sens, list(poolfit(core[core$strict_exact,], "strict_exact"),
                     poolfit(core[core$all_usable,], "all_usable")))
sens <- do.call(rbind, Filter(Negate(is.null), sens))
write.csv(sens, file.path(TAB,"prevalence_regression_sensitivity.csv"), row.names=FALSE)

## ---- 10. Monte Carlo propagation (5000): prevalence-only and joint (prevalence+rg) ----
mc_run <- function(dat, sample_rg, iters=5000) {
  if(nrow(dat)<5) return(NULL)
  bs <- numeric(iters); ps <- numeric(iters)
  for(i in seq_len(iters)) {
    y <- rnorm(nrow(dat), dat$log10prev, dat$se_log10)
    x <- if(sample_rg) rnorm(nrow(dat), dat$rg, ifelse(is.finite(dat$rg_se), dat$rg_se, 0)) else dat$rg
    m <- lm(y ~ x); sm <- summary(m)$coefficients
    bs[i] <- sm["x","Estimate"]; ps[i] <- sm["x","Pr(>|t|)"]
  }
  data.frame(median_beta=median(bs), ci2.5=quantile(bs,.025), ci97.5=quantile(bs,.975),
             prop_beta_gt0=mean(bs>0), prop_p_lt05=mean(ps<0.05), row.names=NULL)
}
mc <- do.call(rbind, lapply(c("Male","Female"), function(sx) {
  dat <- sub("strict_exact", sx)
  rbind(cbind(sex=sx, type="prev_only", mc_run(dat, FALSE)),
        cbind(sex=sx, type="joint_prev_rg", mc_run(dat, TRUE)))
}))
write.csv(mc, file.path(TAB,"prevalence_monte_carlo_summary.csv"), row.names=FALSE)
# NOTE: type=="joint_prev_rg" is a NOISE-INJECTION sensitivity (rg* ~ N(rg,SE^2) added as
# predictor noise); it mechanically attenuates beta and is NOT an EIV correction. The valid
# EIV correction is SIMEX below.

## ---- 10b. SIMEX errors-in-variables correction for rg measurement error (heteroscedastic) ----
# rg has known per-disease SE_rg. SIMEX adds extra noise at levels lambda, fits beta(lambda),
# then extrapolates (quadratic) to lambda = -1 (= zero measurement error). This CORRECTS the
# attenuation rather than adding to it. 95% CI by disease-level nonparametric bootstrap.
simex_beta <- function(dd, lambdas=c(0,0.5,1,1.5,2), B=100){
  dd <- dd[is.finite(dd$rg_se) & dd$rg_se>0, ]
  if(nrow(dd)<5 || sd(dd$rg)==0) return(NA_real_)
  bl <- sapply(lambdas, function(lam) mean(replicate(B, {
    xstar <- dd$rg + sqrt(lam)*rnorm(nrow(dd),0,dd$rg_se)
    coef(lm(dd$log10prev ~ xstar))[2] })))
  ex <- lm(bl ~ lambdas + I(lambdas^2))
  as.numeric(predict(ex, newdata=data.frame(lambdas=-1)))
}
simex <- do.call(rbind, lapply(c("Male","Female"), function(sx){
  dat <- sub("strict_exact", sx)
  base <- simex_beta(dat)
  boot <- replicate(300, { idx <- sample(nrow(dat), replace=TRUE); simex_beta(dat[idx,]) })
  boot <- boot[is.finite(boot)]
  data.frame(sex=sx, n=nrow(dat),
             naive_beta=unname(coef(lm(log10prev~rg,dat))[2]),
             simex_beta=base, ci_lo=quantile(boot,.025,names=FALSE),
             ci_hi=quantile(boot,.975,names=FALSE),
             prop_boot_gt0=mean(boot>0), row.names=NULL)
}))
write.csv(simex, file.path(TAB,"prevalence_simex_rg_eiv.csv"), row.names=FALSE)

## ---- 10c. standardized beta across all 5 datasets x sex (stable-but-noisier vs created-by-mappings) ----
fit_std <- function(dd, dataset, sx){
  if(nrow(dd)<5 || sd(dd$rg)==0) return(NULL)
  z <- data.frame(y=scale(dd$log10prev)[,1], x=scale(dd$rg)[,1])
  m <- lm(y~x, data=z); V <- vcovHC(m,type="HC3")
  b <- coef(m)["x"]; se <- sqrt(V["x","x"]); tc <- qt(.975, df.residual(m))
  data.frame(dataset=dataset, sex=sx, n=nrow(dd), beta_std=unname(b), se_hc3=unname(se),
             ci_lo=unname(b-tc*se), ci_hi=unname(b+tc*se),
             p=unname(coeftest(m,vcov.=V)["x","Pr(>|t|)"]), row.names=NULL)
}
stdtab <- do.call(rbind, unlist(lapply(DS, function(ds)
  lapply(c("Male","Female"), function(sx) fit_std(sub(ds,sx), ds, sx))), recursive=FALSE))
write.csv(stdtab, file.path(TAB,"prevalence_standardized_beta.csv"), row.names=FALSE)

## ---- 11. influence diagnostics on primary (strict_exact per sex) ----
infl <- do.call(rbind, lapply(c("Male","Female"), function(sx) {
  dat <- sub("strict_exact", sx); if(nrow(dat)<5) return(NULL)
  m <- lm(log10prev ~ rg, data=dat); im <- influence.measures(m)
  b_full <- coef(m)["rg"]
  loo <- sapply(seq_len(nrow(dat)), function(j) coef(lm(log10prev~rg, data=dat[-j,]))["rg"])
  data.frame(sex=sx, analysis_disease_id=dat$analysis_disease_id, GBD_cause_name=dat$GBD_cause_name,
             rg=dat$rg, log10prev=dat$log10prev, category=dat$category,
             loo_beta=loo, loo_delta=loo-unname(b_full),
             cooks_d=cooks.distance(m), leverage=hatvalues(m),
             stud_resid=rstudent(m), dfbeta_rg=im$infmat[,"dfb.rg"], stringsAsFactors=FALSE)
}))
# leave-one-category-out betas appended as summary rows
loco <- do.call(rbind, lapply(c("Male","Female"), function(sx) {
  dat <- sub("strict_exact", sx)
  do.call(rbind, lapply(unique(dat$category), function(cc) {
    dd <- dat[dat$category!=cc,]; if(nrow(dd)<5||length(unique(dd$rg))<3) return(NULL)
    data.frame(sex=sx, analysis_disease_id=paste0("LOCO_drop_",cc), GBD_cause_name=paste0("drop ",cc),
               rg=NA, log10prev=NA, category=cc,
               loo_beta=coef(lm(log10prev~rg,data=dd))["rg"], loo_delta=NA,
               cooks_d=NA, leverage=NA, stud_resid=NA, dfbeta_rg=NA, stringsAsFactors=FALSE)
  }))
}))
write.csv(rbind(infl, loco), file.path(TAB,"prevalence_influence_diagnostics.csv"), row.names=FALSE)

## ---- 12. figures: boxplot/violin (equal weight, individual points) ----
mk_box <- function(ds, fname) {
  png(file.path(FIG, fname), width=1500, height=650, res=150)
  par(mfrow=c(1,2), mar=c(4,4,3,1))
  for(sx in c("Male","Female")) {
    dat <- core[core[[ds]] & core$sex==sx, ]
    if(nrow(dat)<4){ plot.new(); title(paste(sx,"n<4")); next }
    dat$grp <- ifelse(dat$rg<0, "rg<0","rg>0")
    ng <- table(dat$grp)
    wt <- wilcox.test(log10prev~grp, data=dat)
    cols <- ifelse(dat$grp=="rg<0", "#377eb8","#e41a1c")
    boxplot(log10prev~grp, data=dat, outline=FALSE, col="grey92", border="grey40",
            ylab="log10 prevalence (per 100k)", xlab="",
            names=c(sprintf("rg<0 (n=%d)", ng["rg<0"]), sprintf("rg>0 (n=%d)", ng["rg>0"])),
            main=sprintf("%s  Wilcoxon P=%.3g", sx, wt$p.value))
    set.seed(1); xj <- jitter(as.numeric(factor(dat$grp)), amount=.12)
    points(xj, dat$log10prev, pch=ifelse(dat$mapping_type=="exact",19,1), col=cols, cex=1.1)
  }
  mtext(paste0("dataset: ",ds," (equal weight; boxplot descriptive only)"), outer=TRUE, line=-1.2, cex=.8)
  dev.off()
}
for(ds in DS) mk_box(ds, paste0("box_",ds,".png"))

## regression scatter (primary strict_exact) with HC3 fit
png(file.path(FIG,"fig_primary_rg_prevalence_strict_exact.png"), width=1500, height=650, res=150)
par(mfrow=c(1,2), mar=c(4,4,3,1))
for(sx in c("Male","Female")) {
  dat <- sub("strict_exact", sx); r <- primary[primary$sex==sx & primary$model=="unadjusted_HC3",]
  plot(dat$rg, dat$log10prev, pch=19, col="#444444",
       xlab=sprintf("disease-fertility rg (%s)", ifelse(sx=="Male","father","mother")),
       ylab="log10 prevalence (per 100k)",
       main=sprintf("%s  n=%d  beta=%.2f  P=%.3g", sx, r$n, r$beta, r$p))
  abline(lm(log10prev~rg, data=dat), col="#e41a1c", lwd=2)
}
dev.off()

cat("DONE. datasets rows=", nrow(core),
    " strict_exact diseases M/F=", sum(core$strict_exact & core$sex=="Male"),
    "/", sum(core$strict_exact & core$sex=="Female"), "\n", sep="")
print(primary[,c("model","sex","n","beta","se_hc3","p","pearson_r","spearman_r")])
cat("\n-- sensitivity (unadjusted only) --\n")
print(sens[sens$model=="unadjusted_HC3", c("dataset","sex","n","beta","p")])
cat("\n-- monte carlo --\n"); print(mc)
