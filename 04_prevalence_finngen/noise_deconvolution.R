#!/usr/bin/env Rscript
## Correction of the per-endpoint antagonistic share for false positives.
## Observed share s = phi*0.5 + (1-phi)*t, with phi = min(0.05*F/n, 1); every gene is
## treated as null, so phi is an upper bound. phi never reaches 1 here (max 0.64/0.59).
## Cut from the manuscript 2026-08-15; kept for the review response.
suppressMessages(library(data.table))

CLASSES <- Sys.getenv("ANTAG_CLASSES",
  "../data/antag50_classes.tsv")
PREV    <- Sys.getenv("SECTION_TABLE",
  "../data/section_table.tsv")
NPERM   <- 10000L
F_GENES <- c(m = 3883, f = 3649)   # genes passing the nominal fertility threshold

a  <- fread(CLASSES)[mode == "nominal"]
st <- fread(PREV, select = c("phenocode", "prev"))
d  <- merge(a, st, by = "phenocode")[prev > 0][, lp := log10(prev)]
stopifnot(nrow(d) == 2754)

set.seed(1)
for (sx in c("m", "f")) {
  n   <- d[[paste0("n_", sx)]]
  s   <- (d[[paste0("pp_", sx)]] + d[[paste0("mm_", sx)]]) / n
  phi <- pmin(0.05 * F_GENES[[sx]] / n, 1)
  stopifnot(all(phi < 1))                      # t undefined at phi = 1
  t   <- (s - 0.5 * phi) / (1 - phi)
  w   <- n * (1 - phi)

  raw <- cor.test(s, d$lp)
  dec <- cor.test(t, d$lp)
  wr  <- cov.wt(cbind(t, d$lp), wt = w, cor = TRUE)$cor[1, 2]
  perm <- replicate(NPERM,
            cov.wt(cbind(t, sample(d$lp)), wt = w, cor = TRUE)$cor[1, 2])
  pperm <- (1 + sum(abs(perm) >= abs(wr))) / (NPERM + 1)

  cat(sprintf(
    "%s | phi median %.3f max %.3f | raw r=%+.3f | deconvolved r=%+.3f P=%.2g | weighted r=%+.3f perm P=%.1g\n",
    c(m = "male", f = "female")[[sx]], median(phi), max(phi),
    raw$estimate, dec$estimate, dec$p.value, wr, pperm))
}
