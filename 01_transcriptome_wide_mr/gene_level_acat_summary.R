#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

mr_dir <- "${MR}"
out_dir <- "${ROOT}"
component_p_threshold <- 0.05
excluded_files <- c(
  "colorectalcancerwu__male.txt.gz",
  "colorectalcancerwu__female.txt.gz",
  "breastcancerwu__female.txt.gz",
  "lungcancerwu__male.txt.gz",
  "lungcancerwu__female.txt.gz",
  "prostatecancerwu__male.txt.gz",
  "GCST90428116Postmenopausalbreastcancer__female.txt.gz"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(mr_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
files <- files[!basename(files) %in% excluded_files]
if (length(files) == 0) {
  stop("No .txt.gz files found in: ", mr_dir)
}

clip_p <- function(p) {
  p <- as.numeric(p)
  p[is.na(p)] <- NA_real_
  p[p <= 0] <- 1e-300
  p[p >= 1] <- 1 - 1e-16
  p
}

acat_from_sum <- function(tan_sum, n) {
  out <- rep(NA_real_, length(tan_sum))
  ok <- !is.na(tan_sum) & !is.na(n) & n > 0
  stat <- tan_sum[ok] / n[ok]
  out[ok] <- 0.5 - atan(stat) / pi
  out[ok & is.infinite(stat) & stat > 0] <- 0
  out
}

consistency_label <- function(prop_same, n_sig) {
  fifelse(n_sig == 0 | is.na(prop_same), "no_sig_components",
    fifelse(prop_same >= 0.70, "consistent",
      fifelse(prop_same <= 0.30, "opposite_dominant", "mixed")))
}

summarize_file <- function(file, effect) {
  cols <- c("tissue", "disease", "fitness", "gene")
  if (effect == "b1") {
    cols <- c(cols, "b1", "se1", "p1")
    beta_col <- "b1"
    se_col <- "se1"
    p_col <- "p1"
    key_cols <- c("gene", "disease")
    context_col <- "fitness"
  } else {
    cols <- c(cols, "b2", "se2", "p2")
    beta_col <- "b2"
    se_col <- "se2"
    p_col <- "p2"
    key_cols <- c("gene", "fitness")
    context_col <- "disease"
  }

  dt <- fread(file, select = cols, showProgress = FALSE)
  setnames(dt, c(beta_col, se_col, p_col), c("beta", "se", "p"))
  dt <- dt[!is.na(gene) & !is.na(p) & !is.na(beta) & !is.na(se)]
  dt[, p_clip := clip_p(p)]
  dt[, acat_tan_sum := tan((0.5 - p_clip) * pi)]
  dt[, valid_meta := is.finite(se) & se > 0 & is.finite(beta)]
  dt[valid_meta == TRUE, w := 1 / (se * se)]
  dt[is.na(w) | !is.finite(w), w := 0]
  dt[, wb := w * beta]
  dt[, wb2 := w * beta * beta]
  dt[, sig_positive := p < component_p_threshold & beta > 0]
  dt[, sig_negative := p < component_p_threshold & beta < 0]

  setorderv(dt, c(key_cols, "p"))
  top <- dt[, .SD[1], by = key_cols]
  top <- top[, c(key_cols, "p", "tissue", context_col, "beta", "se"), with = FALSE]
  setnames(top,
    c("p", "tissue", context_col, "beta", "se"),
    c("min_p", "top_tissue", "top_context", "top_beta", "top_se"))

  agg <- dt[, .(
    n_tests_combined = .N,
    acat_tan_sum = sum(acat_tan_sum, na.rm = TRUE),
    meta_sum_w = sum(w, na.rm = TRUE),
    meta_sum_wb = sum(wb, na.rm = TRUE),
    meta_sum_wb2 = sum(wb2, na.rm = TRUE),
    n_sig_positive = sum(sig_positive, na.rm = TRUE),
    n_sig_negative = sum(sig_negative, na.rm = TRUE)
  ), by = key_cols]

  out <- merge(agg, top, by = key_cols, all.x = TRUE, sort = FALSE)
  out
}

merge_accumulator <- function(acc, part, key_cols) {
  if (is.null(acc)) {
    return(part)
  }

  merged <- merge(acc, part, by = key_cols, all = TRUE, suffixes = c("_old", "_new"), sort = FALSE)

  sum_cols <- c(
    "n_tests_combined", "acat_tan_sum", "meta_sum_w", "meta_sum_wb",
    "meta_sum_wb2", "n_sig_positive", "n_sig_negative"
  )
  for (col in sum_cols) {
    old <- paste0(col, "_old")
    new <- paste0(col, "_new")
    merged[, (col) := fcoalesce(as.numeric(get(old)), 0) + fcoalesce(as.numeric(get(new)), 0)]
  }

  use_new <- is.na(merged$min_p_old) |
    (!is.na(merged$min_p_new) & merged$min_p_new < merged$min_p_old)

  for (col in c("min_p", "top_tissue", "top_context", "top_beta", "top_se")) {
    old <- paste0(col, "_old")
    new <- paste0(col, "_new")
    merged[, (col) := fifelse(use_new, get(new), get(old))]
  }

  keep <- c(key_cols, sum_cols, "min_p", "top_tissue", "top_context", "top_beta", "top_se")
  merged[, keep, with = FALSE]
}

finalize_summary <- function(acc, key_cols, target_col, effect_name) {
  acc[, p_acat := acat_from_sum(acat_tan_sum, n_tests_combined)]
  acc[, p_bonf_across_genes := p.adjust(p_acat, method = "bonferroni"), by = target_col]
  acc[, p_bonf_global := p.adjust(p_acat, method = "bonferroni")]
  acc[, primary_direction := fifelse(top_beta > 0, "positive",
    fifelse(top_beta < 0, "negative", "zero"))]

  acc[, n_sig_components := n_sig_positive + n_sig_negative]
  acc[, n_same_sign_sig := fifelse(primary_direction == "positive", n_sig_positive,
    fifelse(primary_direction == "negative", n_sig_negative, 0))]
  acc[, n_opposite_sign_sig := fifelse(primary_direction == "positive", n_sig_negative,
    fifelse(primary_direction == "negative", n_sig_positive, 0))]
  acc[, sign_prop_same_sig := fifelse(n_sig_components > 0,
    n_same_sign_sig / n_sig_components, NA_real_)]
  acc[, sign_consistency := consistency_label(sign_prop_same_sig, n_sig_components)]

  acc[, meta_beta_ivw := fifelse(meta_sum_w > 0, meta_sum_wb / meta_sum_w, NA_real_)]
  acc[, meta_se_ivw := fifelse(meta_sum_w > 0, sqrt(1 / meta_sum_w), NA_real_)]
  acc[, meta_z_ivw := meta_beta_ivw / meta_se_ivw]
  acc[, meta_p_ivw := 2 * pnorm(-abs(meta_z_ivw))]

  acc[, heterogeneity_q := fifelse(
    meta_sum_w > 0,
    pmax(0, meta_sum_wb2 - (meta_sum_wb * meta_sum_wb / meta_sum_w)),
    NA_real_
  )]
  acc[, heterogeneity_df := pmax(n_tests_combined - 1, 0)]
  acc[, heterogeneity_p := fifelse(
    heterogeneity_df > 0,
    pchisq(heterogeneity_q, df = heterogeneity_df, lower.tail = FALSE),
    NA_real_
  )]
  acc[, heterogeneity_i2 := fifelse(
    heterogeneity_q > 0 & heterogeneity_df > 0,
    pmax(0, (heterogeneity_q - heterogeneity_df) / heterogeneity_q),
    0
  )]

  acc[, effect := effect_name]
  acc[, component_p_threshold := component_p_threshold]

  out_cols <- c(
    "gene", target_col, "effect", "p_acat", "p_bonf_across_genes",
    "p_bonf_global", "n_tests_combined",
    "min_p", "top_tissue", "top_context", "top_beta", "top_se",
    "primary_direction", "n_sig_components", "n_same_sign_sig",
    "n_opposite_sign_sig", "sign_prop_same_sig", "sign_consistency",
    "meta_beta_ivw", "meta_se_ivw", "meta_z_ivw", "meta_p_ivw",
    "heterogeneity_q", "heterogeneity_df", "heterogeneity_p",
    "heterogeneity_i2", "component_p_threshold"
  )
  setorderv(acc, c(target_col, "p_acat", "gene"))
  acc[, out_cols, with = FALSE]
}

trait_summary <- function(b1_out, b2_out) {
  b1_summary <- b1_out[, .(
    effect = "b1_disease",
    trait = disease[1],
    n_genes_tested = uniqueN(gene),
    n_sig_bonf_across_genes = sum(p_bonf_across_genes < 0.05, na.rm = TRUE),
    n_sig_bonf_global = sum(p_bonf_global < 0.05, na.rm = TRUE),
    n_sig_consistent_across_genes = sum(
      p_bonf_across_genes < 0.05 & sign_consistency == "consistent",
      na.rm = TRUE
    ),
    n_sig_mixed_across_genes = sum(
      p_bonf_across_genes < 0.05 & sign_consistency == "mixed",
      na.rm = TRUE
    ),
    n_sig_positive_across_genes = sum(
      p_bonf_across_genes < 0.05 & primary_direction == "positive",
      na.rm = TRUE
    ),
    n_sig_negative_across_genes = sum(
      p_bonf_across_genes < 0.05 & primary_direction == "negative",
      na.rm = TRUE
    )
  ), by = disease][, disease := NULL]

  b2_summary <- b2_out[, .(
    effect = "b2_fitness",
    trait = fitness[1],
    n_genes_tested = uniqueN(gene),
    n_sig_bonf_across_genes = sum(p_bonf_across_genes < 0.05, na.rm = TRUE),
    n_sig_bonf_global = sum(p_bonf_global < 0.05, na.rm = TRUE),
    n_sig_consistent_across_genes = sum(
      p_bonf_across_genes < 0.05 & sign_consistency == "consistent",
      na.rm = TRUE
    ),
    n_sig_mixed_across_genes = sum(
      p_bonf_across_genes < 0.05 & sign_consistency == "mixed",
      na.rm = TRUE
    ),
    n_sig_positive_across_genes = sum(
      p_bonf_across_genes < 0.05 & primary_direction == "positive",
      na.rm = TRUE
    ),
    n_sig_negative_across_genes = sum(
      p_bonf_across_genes < 0.05 & primary_direction == "negative",
      na.rm = TRUE
    )
  ), by = fitness][, fitness := NULL]

  out <- rbind(b1_summary, b2_summary, use.names = TRUE)
  setorder(out, effect, -n_sig_bonf_across_genes, trait)
  out
}

message("Input files after exclusions: ", length(files))
message("Excluded files: ", paste(excluded_files, collapse = ", "))

b1_acc <- NULL
b2_acc <- NULL

for (i in seq_along(files)) {
  file <- files[i]
  message(sprintf("[%03d/%03d] %s", i, length(files), basename(file)))

  b1_part <- summarize_file(file, "b1")
  b1_acc <- merge_accumulator(b1_acc, b1_part, c("gene", "disease"))
  rm(b1_part)
  gc(FALSE)

  b2_part <- summarize_file(file, "b2")
  b2_acc <- merge_accumulator(b2_acc, b2_part, c("gene", "fitness"))
  rm(b2_part)
  gc(FALSE)
}

b1_out <- finalize_summary(b1_acc, c("gene", "disease"), "disease", "b1_disease")
b2_out <- finalize_summary(b2_acc, c("gene", "fitness"), "fitness", "b2_fitness")

b1_file <- file.path(out_dir, "gene_level_b1_disease_acat_summary.csv")
b2_file <- file.path(out_dir, "gene_level_b2_fitness_acat_summary.csv")

fwrite(b1_out, b1_file)
fwrite(b2_out, b2_file)
fwrite(
  trait_summary(b1_out, b2_out),
  file.path(out_dir, "significant_genes_per_trait_summary.csv")
)

manifest <- data.table(
  output = c(
    "gene_level_b1_disease_acat_summary.csv",
    "gene_level_b2_fitness_acat_summary.csv",
    "significant_genes_per_trait_summary.csv"
  ),
  rows = c(nrow(b1_out), nrow(b2_out), nrow(trait_summary(b1_out, b2_out))),
  description = c(
    "Gene-disease summary of b1 disease-effect p1 aggregated by ACAT across tissues and fitness traits.",
    "Gene-fitness summary of b2 fitness-effect p2 aggregated by ACAT across tissues and diseases.",
    "Per-trait counts of Bonferroni-significant genes from the gene-level ACAT summaries."
  ),
  component_p_threshold = component_p_threshold,
  bonferroni_across_genes_scope = c(
    "Bonferroni-adjusted p_acat across genes within each disease",
    "Bonferroni-adjusted p_acat across genes within each fitness trait",
    "Counts use p_bonf_across_genes < 0.05"
  ),
  bonferroni_global_scope = c(
    "Bonferroni-adjusted p_acat over all gene x disease tests",
    "Bonferroni-adjusted p_acat over all gene x fitness tests",
    "Counts also report p_bonf_global < 0.05"
  ),
  excluded_files = paste(excluded_files, collapse = ";"),
  mr_result_dir = mr_dir
)
fwrite(manifest, file.path(out_dir, "gene_level_acat_summary_manifest.csv"))

message("Wrote: ", b1_file)
message("Wrote: ", b2_file)
message("Done.")
