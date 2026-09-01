#!/usr/bin/env Rscript
## The single pre-specified model behind the "Disease-fertility coupling predicts
## prevalence beyond mortality and disease onset" section. Every partial R^2 in
## that section and in Fig. 4d-f comes from this one fit -- the model is never
## re-specified to make a term look larger.
##
##   log10 prevalence ~ coupling with male fertility
##                    + coupling with female fertility
##                    + log 5-year mortality hazard ratio
##                    + ns(median age at first recorded event, 3)
##                    + ICD chapter
##
## Age at first recorded event enters as a natural spline: a linear term averages
## the two opposing arms of a non-monotonic relationship to near zero
## (dR^2 = 0.028, P = 1e-17 for spline over linear).
##
## Inputs
##   model_input_1678.tsv   phenocode, lp (= log10 prevalence %), chap (ICD chapter,
##                          20 levels), rho_GE_male, rho_GE_female,
##                          lmort (= log mean of female and male 5-year mortality HR),
##                          onset_all (Risteys median age at first recorded event)
## Sources: rho from fusios_full_20260824/tables_fixed/rhoGE_50tissue_matched.tsv;
##          prevalence, onset and mortality from Risteys R12 (risteys_r12_epi.tsv).
suppressMessages({library(data.table); library(splines)})
d <- fread("model_input_1678.tsv")

FULL <- lp ~ factor(chap) + rho_GE_male + rho_GE_female + lmort + ns(onset_all, 3)
mf   <- lm(FULL, d)
cat(sprintf("n = %d   coefficients = %d   residual df = %d   R2 = %.4f   adj R2 = %.4f\n\n",
            nrow(d), length(coef(mf)), df.residual(mf),
            summary(mf)$r.squared, summary(mf)$adj.r.squared))

## partial R^2 = (R2_full - R2_reduced) / (1 - R2_reduced); cross-checked against
## t^2 / (t^2 + residual df) for the one-degree-of-freedom terms.
pr2 <- function(drop) {
  m0 <- update(mf, as.formula(paste(". ~ . -", drop)))
  c(partial_R2 = (summary(mf)$r.squared - summary(m0)$r.squared) / (1 - summary(m0)$r.squared),
    dR2        = summary(mf)$r.squared - summary(m0)$r.squared,
    P          = anova(m0, mf)[2, "Pr(>F)"])
}
terms <- list(
  "coupling, male fertility"        = "rho_GE_male",
  "coupling, female fertility"      = "rho_GE_female",
  "coupling, joint (2 df)"          = "rho_GE_male - rho_GE_female",
  "log 5-year mortality HR"         = "lmort",
  "age at first recorded event (3 df)" = "ns(onset_all, 3)",
  "ICD chapter (19 df)"             = "factor(chap)")
for (nm in names(terms)) {
  v <- pr2(terms[[nm]])
  cat(sprintf("  %-36s partial R2 = %.4f   dR2 = %.4f   P = %.2g\n",
              nm, v["partial_R2"], v["dR2"], v["P"]))
}

## spline vs linear, inside the same model
lin <- lm(lp ~ factor(chap) + rho_GE_male + rho_GE_female + lmort + onset_all, d)
cat(sprintf("\nspline over linear onset: dR2 = %.4f   P = %.2g\n",
            summary(mf)$r.squared - summary(lin)$r.squared, anova(lin, mf)[2, "Pr(>F)"]))

## model-free description used in Fig. 4f: coupling-prevalence correlation within
## tertiles of age at first recorded event
d[, ot := cut(onset_all, quantile(onset_all, 0:3/3), include.lowest = TRUE,
              labels = c("early", "mid", "late"))]
cat("\ncoupling vs log10 prevalence within onset tertiles:\n")
print(d[!is.na(ot), .(n = .N, median_onset = round(median(onset_all), 1),
                      r_male   = round(cor(rho_GE_male, lp), 3),
                      r_female = round(cor(rho_GE_female, lp), 3)), by = ot][order(ot)])

## Covariates deliberately excluded: case count, genome-wide significant locus count
## and the number of genes contributing to each coupling estimate. All three are
## downstream of the disease case count, which is nearly collinear with prevalence
## (r = 0.97 on the log10 scale), so conditioning on them removes the outcome
## variance the model is asked to explain.
