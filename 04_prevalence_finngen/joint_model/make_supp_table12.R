suppressMessages({library(openxlsx); library(data.table); library(splines)})
setwd("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait")

## ---- build the joint-model coefficient table ----
d <- fread("model_input_1678.tsv")
FULL <- lp ~ factor(chap) + rho_GE_male + rho_GE_female + lmort + ns(onset_all, 3)
mf <- lm(FULL, d)
pr <- function(drop) {
  m0 <- update(mf, as.formula(paste(". ~ . -", drop)))
  data.table(partial_R2 = (summary(mf)$r.squared - summary(m0)$r.squared) / (1 - summary(m0)$r.squared),
             delta_R2   = summary(mf)$r.squared - summary(m0)$r.squared,
             df         = abs(df.residual(mf) - df.residual(m0)),
             P          = anova(m0, mf)[2, "Pr(>F)"])
}
z <- copy(d); for (v in c("rho_GE_male","rho_GE_female","lmort","onset_all")) z[[v]] <- as.numeric(scale(z[[v]]))
z$lp <- as.numeric(scale(z$lp))
bz <- summary(lm(lp ~ factor(chap) + rho_GE_male + rho_GE_female + lmort + onset_all, z))$coefficients
rows <- list(
  c("Coupling with male fertility",        "rho_GE_male"),
  c("Coupling with female fertility",      "rho_GE_female"),
  c("Coupling, both terms jointly",        "rho_GE_male - rho_GE_female"),
  c("log 5-year mortality hazard ratio",   "lmort"),
  c("Median age at first recorded event (natural spline, 3 df)", "ns(onset_all, 3)"),
  c("ICD chapter (19 df)",                 "factor(chap)"))
out <- rbindlist(lapply(rows, function(r) {
  v <- pr(r[2])
  key <- gsub("ns\\(onset_all, 3\\)", "onset_all", r[2])
  beta <- if (key %in% rownames(bz)) bz[key, 1] else NA_real_
  se   <- if (key %in% rownames(bz)) bz[key, 2] else NA_real_
  data.table(Term = r[1], df = v$df,
             Standardised_beta = round(beta, 4), SE = round(se, 4),
             partial_R2 = round(v$partial_R2, 4), delta_R2 = round(v$delta_R2, 4),
             P = signif(v$P, 3))
}))
lin <- lm(lp ~ factor(chap) + rho_GE_male + rho_GE_female + lmort + onset_all, d)
out <- rbind(out, data.table(Term = "Onset spline over a linear onset term", df = 2,
  Standardised_beta = NA_real_, SE = NA_real_,
  partial_R2 = NA_real_,
  delta_R2 = round(summary(mf)$r.squared - summary(lin)$r.squared, 4),
  P = signif(anova(lin, mf)[2, "Pr(>F)"], 3)))

desc <- paste0("Supplementary Table 12. Joint model of population prevalence across the 1,678 FinnGen ",
"endpoints with a Risteys five-year mortality estimate. Model: log10 prevalence ~ coupling with male ",
"fertility + coupling with female fertility + log five-year mortality hazard ratio + natural spline ",
"(3 df) in the median age at first recorded event + ICD chapter (20 levels). n = 1,678; R2 = ",
sprintf("%.3f", summary(mf)$r.squared), "; adjusted R2 = ", sprintf("%.3f", summary(mf)$adj.r.squared),
"; residual df = ", df.residual(mf), ". partial R2 = (R2_full - R2_reduced)/(1 - R2_reduced) from ",
"dropping the term from this one model; delta R2 is the same difference on the total-variance scale. ",
"Standardised betas come from the same model with a linear onset term so that every predictor has a ",
"single coefficient. Case counts, genome-wide significant locus counts and the number of genes ",
"behind each coupling estimate are deliberately excluded: all are downstream of the disease case ",
"count, which is nearly collinear with prevalence (r = 0.97 on the log10 scale).")

wb <- loadWorkbook("Supplementary Table.xlsx")
if ("Supplementary Table 12" %in% names(wb)) removeWorksheet(wb, "Supplementary Table 12")
addWorksheet(wb, "Supplementary Table 12")
writeData(wb, "Supplementary Table 12", desc, startRow = 1, colNames = FALSE)
writeData(wb, "Supplementary Table 12", out, startRow = 2)

## restore numeric sheet order (Contents first)
nm  <- names(wb); num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", nm)))
num[nm == "Contents"] <- -1
worksheetOrder(wb) <- order(num)
saveWorkbook(wb, "Supplementary Table.xlsx", overwrite = TRUE)
print(out)
cat("\n工作表顺序:", paste(names(wb)[worksheetOrder(wb)], collapse = " | "), "\n")
