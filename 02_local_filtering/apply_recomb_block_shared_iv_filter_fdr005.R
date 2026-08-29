suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- "${HOME}/Library/CloudStorage/Dropbox/Lin/Jianhai/conflict_genes/noffspring_rerun/mr_res_97_acat_summary"
work_dir <- "${ROOT}"
shared_dir <- file.path(work_dir, "results/shared_iv_cluster_filter")
dir.create(file.path(work_dir, "data"), recursive = TRUE, showWarnings = FALSE)

b1_file <- file.path(base_dir, "gene_level_b1_disease_acat_summary.csv")
b2_file <- file.path(base_dir, "gene_level_b2_fitness_acat_summary.csv")
pair_file <- file.path(shared_dir, "fdr005_gene_level_shared_iv_pair_metrics.csv.gz")
gtf_file <- file.path(work_dir, "data/gencode.v19.annotation.gtf.gz")
gene_coords_file <- file.path(work_dir, "data/gene_coords_v19.tsv")
recomb_dir <- "${DATA}/plink.GRCh37.map"

significance_method <- Sys.getenv("SIGNIFICANCE_METHOD", "fdr")
sig_cutoff <- as.numeric(Sys.getenv("SIG_CUTOFF", Sys.getenv("FDR_CUTOFF", "0.05")))
mixed_artifact_mode <- Sys.getenv("MIXED_ARTIFACT_MODE", "confirmed")
sig_tag <- if (significance_method == "bonf_global") {
  sprintf("bonf%03d", as.integer(round(sig_cutoff * 1000)))
} else {
  sprintf("fdr%03d", as.integer(round(sig_cutoff * 1000)))
}
if (mixed_artifact_mode == "broad") sig_tag <- paste0(sig_tag, "_broadmixed")
out_dir <- file.path(work_dir, sprintf("results/recomb_block_shared_iv_filter_%s", sig_tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cm_gap_max <- 0.5
shared_n_min <- 5L
shared_fraction_min <- 0.5
same_sign_eqtl_cor_min <- 0.8
flipper_eqtl_cor_max <- -0.8
mhc_chr <- 6L
mhc_start <- 25000000L
mhc_end <- 34000000L

required <- c(b1_file, b2_file, pair_file, gtf_file, recomb_dir)
missing_required <- required[!file.exists(required)]
if (length(missing_required) > 0) {
  stop("Missing required input(s):\n", paste(missing_required, collapse = "\n"))
}

extract_attr <- function(x, key) {
  pat <- paste0(key, " \"([^\"]+)\"")
  hit <- regmatches(x, regexpr(pat, x))
  sub(pat, "\\1", hit)
}

find_components <- function(nodes, edges) {
  nodes <- sort(unique(nodes))
  if (length(nodes) == 0) return(data.table(gene = character(), component = integer()))
  parent <- setNames(nodes, nodes)
  find_root <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[[rb]] <<- ra
  }
  if (nrow(edges) > 0) {
    for (i in seq_len(nrow(edges))) union_nodes(edges$gene_a[i], edges$gene_b[i])
  }
  data.table(gene = nodes, root = vapply(nodes, find_root, character(1)))[
    , component := .GRP, by = root
  ][, .(gene, component)]
}

if (!file.exists(gene_coords_file)) {
  message("Creating gene coordinate file from GENCODE v19")
  gtf <- fread(
    gtf_file,
    sep = "\t",
    header = FALSE,
    quote = "",
    col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attrs")
  )
  genes <- gtf[feature == "gene" & grepl("^chr([0-9]+|X|Y)$", chr)]
  genes[, gene := sub("\\..*$", "", extract_attr(attrs, "gene_id"))]
  genes[, chr_num := suppressWarnings(as.integer(sub("^chr", "", chr)))]
  genes <- genes[!is.na(chr_num) & chr_num >= 1 & chr_num <= 22]
  genes[, pos := fifelse(strand == "-", end, start)]
  coords <- unique(genes[, .(gene, chr = chr_num, pos, start, end, strand)])
  fwrite(coords[order(chr, pos, gene)], gene_coords_file, sep = "\t")
}

message("Reading recombination maps")
recomb <- rbindlist(lapply(1:22, function(chr_i) {
  f <- file.path(recomb_dir, sprintf("plink.chr%d.GRCh37.map", chr_i))
  x <- fread(f, header = FALSE, col.names = c("chr", "id", "cm", "bp"))
  x[, chr := chr_i]
  x[, .(chr, bp, cm)]
}), use.names = TRUE)
setorder(recomb, chr, bp)

coords <- fread(gene_coords_file)
coords[, chr := as.integer(chr)]
coords[, pos := as.integer(pos)]
coords <- coords[chr %between% c(1L, 22L)]
coords[, cm := NA_real_]
for (chr_i in 1:22) {
  map_i <- recomb[chr == chr_i]
  idx <- which(coords$chr == chr_i)
  if (length(idx) > 0) {
    coords$cm[idx] <- approx(map_i$bp, map_i$cm, xout = coords$pos[idx], rule = 2, ties = "ordered")$y
  }
}
coords[, in_mhc := chr == mhc_chr & pos >= mhc_start & pos <= mhc_end]
fwrite(coords[order(chr, pos, gene)], file.path(out_dir, "gene_coords_v19_with_cm.csv"))

message("Reading MR summaries and building combo universe")
b1 <- fread(
  b1_file,
  select = c("gene", "disease", "p_acat", "p_bonf_global", "meta_beta_ivw", "meta_se_ivw", "meta_z_ivw", "top_beta", "top_se")
)
b1[, b1_q_fdr_check := p.adjust(p_acat, method = "BH")]
b1[, disease_sign := fifelse(meta_beta_ivw > 0, "+", fifelse(meta_beta_ivw < 0, "-", NA_character_))]
b1[, disease_strength := fifelse(!is.na(meta_z_ivw), abs(meta_z_ivw), abs(top_beta / top_se))]

b2 <- fread(
  b2_file,
  select = c("gene", "fitness", "p_acat", "p_bonf_global", "meta_beta_ivw", "meta_z_ivw")
)
b2[, b2_q_fdr_check := p.adjust(p_acat, method = "BH")]

if (significance_method == "bonf_global") {
  b1_sig <- b1[p_bonf_global < sig_cutoff]
  b2_sig <- b2[p_bonf_global < sig_cutoff]
} else if (significance_method == "fdr") {
  b1_sig <- b1[b1_q_fdr_check < sig_cutoff]
  b2_sig <- b2[b2_q_fdr_check < sig_cutoff]
} else {
  stop("Unsupported SIGNIFICANCE_METHOD: ", significance_method)
}

b1_sig <- b1_sig[, .(
  gene,
  disease,
  b1 = meta_beta_ivw,
  b1_q_fdr = b1_q_fdr_check
)]
b2_sig <- b2_sig[, .(
  gene,
  fitness,
  b2 = meta_beta_ivw,
  b2_q_fdr = b2_q_fdr_check
)]
combo <- merge(b1_sig, b2_sig, by = "gene", allow.cartesian = TRUE)
combo[, combo := paste0(
  fifelse(b1 > 0, "+", fifelse(b1 < 0, "-", "0")),
  fifelse(b2 > 0, "+", fifelse(b2 < 0, "-", "0"))
)]
combo <- combo[combo %chin% c("++", "+-", "-+", "--")]
analysis_genes <- sort(unique(combo$gene))
b1 <- b1[gene %in% analysis_genes]
b2 <- b2[gene %in% analysis_genes]

combo <- merge(
  combo,
  b1[, .(gene, disease, disease_sign, disease_strength, b1_meta_beta = meta_beta_ivw, b1_meta_z = meta_z_ivw)],
  by = c("gene", "disease"),
  all.x = TRUE
)
combo <- merge(
  combo,
  coords[, .(gene, chr, pos, cm, in_mhc)],
  by = "gene",
  all.x = TRUE
)
combo[, missing_coords := is.na(chr) | is.na(pos) | is.na(cm)]

message("Reading pair metrics")
pairs <- fread(pair_file)
pairs[, pair_key := fifelse(gene_a < gene_b, paste(gene_a, gene_b, sep = "||"), paste(gene_b, gene_a, sep = "||"))]
pairs[, high_shared_iv := shared_n >= shared_n_min & shared_fraction_min_iv >= shared_fraction_min]
pairs[, redundant_same_sign_pair := high_shared_iv == TRUE & !is.na(eqtl_cor) & eqtl_cor >= same_sign_eqtl_cor_min]
pairs[, artifact_flipper_pair := high_shared_iv == TRUE & !is.na(eqtl_cor) & eqtl_cor <= flipper_eqtl_cor_max]
pair_lookup <- pairs[, .(
  pair_key,
  gene_a,
  gene_b,
  shared_n,
  n_iv_a,
  n_iv_b,
  shared_fraction_min_iv,
  eqtl_cor,
  high_shared_iv,
  redundant_same_sign_pair,
  artifact_flipper_pair
)]

message("Applying recombination-block plus shared-IV/eQTL filter")
combo[, `:=`(
  recomb_block_id = NA_character_,
  recomb_block_n_genes = NA_integer_,
  recomb_block_disease_sign_pattern = NA_character_,
  filter_action = NA_character_,
  filter_reason = NA_character_,
  representative_gene = NA_character_
)]

block_summaries <- list()
pair_evidence_rows <- list()
annotated_parts <- vector("list", 0)
part_i <- 0L
block_i <- 0L
pair_i <- 0L

setorder(combo, disease, fitness, chr, cm, pos, gene)
strata <- unique(combo[, .(disease, fitness)])

for (s in seq_len(nrow(strata))) {
  disease_s <- strata$disease[s]
  fitness_s <- strata$fitness[s]
  dt <- combo[disease == disease_s & fitness == fitness_s]

  no_coord <- dt[missing_coords == TRUE]
  if (nrow(no_coord) > 0) {
    no_coord[, `:=`(
      filter_action = "keep_no_gene_coord",
      filter_reason = "No GRCh37 gene coordinate/cM available; retained because block-level redundancy cannot be evaluated",
      representative_gene = gene
    )]
    part_i <- part_i + 1L
    annotated_parts[[part_i]] <- no_coord
  }

  mhc <- dt[missing_coords == FALSE & in_mhc == TRUE]
  if (nrow(mhc) > 0) {
    mhc[, `:=`(
      recomb_block_id = paste(disease_s, fitness_s, "MHC", sep = "||"),
      recomb_block_n_genes = .N,
      recomb_block_disease_sign_pattern = "MHC",
      filter_action = "drop_mhc",
      filter_reason = "MHC chr6:25-34Mb excluded before recombination-block filtering",
      representative_gene = NA_character_
    )]
    part_i <- part_i + 1L
    annotated_parts[[part_i]] <- mhc
  }

  x <- dt[missing_coords == FALSE & in_mhc == FALSE]
  if (nrow(x) == 0) next
  setorder(x, chr, cm, pos, gene)
  x[, block_index := cumsum(is.na(shift(chr)) | chr != shift(chr) | (cm - shift(cm)) > cm_gap_max)]

  for (b in sort(unique(x$block_index))) {
    block_i <- block_i + 1L
    xb <- copy(x[block_index == b])
    block_id <- paste(disease_s, fitness_s, sprintf("block%05d", block_i), sep = "||")
    genes_b <- sort(unique(xb$gene))
    signs <- sort(unique(na.omit(xb$disease_sign)))
    sign_pattern <- if (length(signs) == 1) {
      paste0("same_", signs)
    } else if (length(signs) > 1) {
      "mixed"
    } else {
      "missing_sign"
    }
    xb[, `:=`(
      recomb_block_id = block_id,
      recomb_block_n_genes = length(genes_b),
      recomb_block_disease_sign_pattern = sign_pattern,
      filter_action = "keep",
      filter_reason = "No block-level redundancy/artifact evidence",
      representative_gene = gene
    )]

    block_pair <- data.table()
    if (length(genes_b) >= 2) {
      gene_pairs <- CJ(gene_a = genes_b, gene_b = genes_b)[gene_a < gene_b]
      gene_pairs[, pair_key := paste(gene_a, gene_b, sep = "||")]
      block_pair <- merge(gene_pairs, pair_lookup, by = c("pair_key", "gene_a", "gene_b"), all.x = TRUE)
      block_pair[is.na(high_shared_iv), `:=`(
        shared_n = 0L,
        high_shared_iv = FALSE,
        redundant_same_sign_pair = FALSE,
        artifact_flipper_pair = FALSE
      )]
      signs_pair <- xb[, .(gene, disease_sign)]
      block_pair <- merge(block_pair, signs_pair[, .(gene_a = gene, sign_a = disease_sign)], by = "gene_a", all.x = TRUE)
      block_pair <- merge(block_pair, signs_pair[, .(gene_b = gene, sign_b = disease_sign)], by = "gene_b", all.x = TRUE)
      block_pair[, opposite_disease_sign := !is.na(sign_a) & !is.na(sign_b) & sign_a != sign_b]
      block_pair[, same_disease_sign := !is.na(sign_a) & !is.na(sign_b) & sign_a == sign_b]
      pair_i <- pair_i + 1L
      pair_evidence_rows[[pair_i]] <- block_pair[, .(
        disease = disease_s,
        fitness = fitness_s,
        recomb_block_id = block_id,
        gene_a,
        gene_b,
        sign_a,
        sign_b,
        shared_n,
        n_iv_a,
        n_iv_b,
        shared_fraction_min_iv,
        eqtl_cor,
        high_shared_iv,
        redundant_same_sign_pair,
        artifact_flipper_pair,
        same_disease_sign,
        opposite_disease_sign
      )]
    }

    if (length(genes_b) == 1) {
      xb[, `:=`(
        filter_action = "keep_singleton_block",
        filter_reason = "Singleton recombination block"
      )]
    } else if (sign_pattern != "mixed") {
      red_edges <- block_pair[redundant_same_sign_pair == TRUE & same_disease_sign == TRUE, .(gene_a, gene_b)]
      if (nrow(red_edges) > 0) {
        comps <- find_components(genes_b, red_edges)
        xb <- merge(xb, comps, by = "gene", all.x = TRUE)
        redundant_components <- comps[, .N, by = component][N > 1, component]
        for (cc in redundant_components) {
          members <- xb[component == cc]
          rep_gene <- members[order(-disease_strength, b1_q_fdr, gene)]$gene[1]
          xb[component == cc, representative_gene := rep_gene]
          xb[component == cc & gene == rep_gene, `:=`(
            filter_action = "keep_representative_same_sign_redundant",
            filter_reason = sprintf(
              "Same-sign block; representative retained for shared-IV redundant subgroup (shared_n>=%d, shared_fraction>=%.2f, eqtl_cor>=%.2f)",
              shared_n_min, shared_fraction_min, same_sign_eqtl_cor_min
            )
          )]
          xb[component == cc & gene != rep_gene, `:=`(
            filter_action = "drop_redundant_same_sign",
            filter_reason = sprintf(
              "Same-sign block; redundant with representative by shared IVs and positive eQTL correlation (shared_n>=%d, shared_fraction>=%.2f, eqtl_cor>=%.2f)",
              shared_n_min, shared_fraction_min, same_sign_eqtl_cor_min
            )
          )]
        }
        xb[is.na(filter_action) | filter_action == "keep", `:=`(
          filter_action = "keep_same_sign_not_redundant",
          filter_reason = "Same-sign block but no supported shared-IV/eQTL redundancy"
        )]
        xb[, component := NULL]
      } else {
        xb[, `:=`(
          filter_action = "keep_same_sign_not_redundant",
          filter_reason = "Same-sign block but no supported shared-IV/eQTL redundancy"
        )]
      }
    } else {
      confirmed_edges <- block_pair[artifact_flipper_pair == TRUE & opposite_disease_sign == TRUE, .(gene_a, gene_b)]
      possible_edges <- block_pair[high_shared_iv == TRUE & opposite_disease_sign == TRUE & artifact_flipper_pair != TRUE, .(gene_a, gene_b)]
      drop_edges <- if (mixed_artifact_mode == "broad") {
        rbindlist(list(
          confirmed_edges[, artifact_label := "confirmed_mixed_artifact"],
          possible_edges[, artifact_label := "possible_mixed_artifact"]
        ), use.names = TRUE, fill = TRUE)
      } else {
        confirmed_edges[, artifact_label := "confirmed_mixed_artifact"]
      }
      drop_edges <- unique(drop_edges)

      if (nrow(drop_edges) > 0) {
        comps <- find_components(genes_b, drop_edges[, .(gene_a, gene_b)])
        xb <- merge(xb, comps, by = "gene", all.x = TRUE)
        edge_components <- merge(drop_edges, comps[, .(gene_a = gene, component)], by = "gene_a", all.x = TRUE)
        component_labels <- edge_components[
          ,
          .(
            has_confirmed = any(artifact_label == "confirmed_mixed_artifact"),
            has_possible = any(artifact_label == "possible_mixed_artifact")
          ),
          by = component
        ]
        confirmed_components <- component_labels[has_confirmed == TRUE, component]
        possible_components <- component_labels[has_confirmed == FALSE & has_possible == TRUE, component]

        xb[component %in% confirmed_components, `:=`(
          filter_action = "drop_confirmed_mixed_artifact",
          filter_reason = sprintf(
            "Mixed-sign block; opposite-sign subgroup shares IVs and has anti-correlated eQTL effects (shared_n>=%d, shared_fraction>=%.2f, eqtl_cor<=%.2f)",
            shared_n_min, shared_fraction_min, flipper_eqtl_cor_max
          ),
          representative_gene = NA_character_
        )]
        xb[component %in% possible_components, `:=`(
          filter_action = "drop_possible_mixed_artifact",
          filter_reason = sprintf(
            "Mixed-sign block; opposite-sign subgroup shares IVs but lacks confirmed anti-correlated eQTL evidence (shared_n>=%d, shared_fraction>=%.2f)",
            shared_n_min, shared_fraction_min
          ),
          representative_gene = NA_character_
        )]
        xb[!(component %in% c(confirmed_components, possible_components)), `:=`(
          filter_action = "keep_mixed_sign_no_shared_iv",
          filter_reason = "Mixed-sign block but no high-shared-IV opposite-sign subgroup"
        )]
        xb[, component := NULL]
      } else {
        xb[, `:=`(
          filter_action = "keep_mixed_sign_no_shared_iv",
          filter_reason = "Mixed-sign block but no high-shared-IV opposite-sign subgroup"
        )]
      }
    }

    block_summaries[[block_i]] <- data.table(
      disease = disease_s,
      fitness = fitness_s,
      recomb_block_id = block_id,
      chr = xb$chr[1],
      min_pos = min(xb$pos, na.rm = TRUE),
      max_pos = max(xb$pos, na.rm = TRUE),
      min_cm = min(xb$cm, na.rm = TRUE),
      max_cm = max(xb$cm, na.rm = TRUE),
      n_genes = length(genes_b),
      disease_sign_pattern = sign_pattern,
      n_redundant_same_sign_pairs = if (nrow(block_pair)) sum(block_pair$redundant_same_sign_pair == TRUE & block_pair$same_disease_sign == TRUE, na.rm = TRUE) else 0L,
      n_artifact_flipper_pairs = if (nrow(block_pair)) sum(block_pair$artifact_flipper_pair == TRUE & block_pair$opposite_disease_sign == TRUE, na.rm = TRUE) else 0L,
      n_drop_redundant = sum(xb$filter_action == "drop_redundant_same_sign", na.rm = TRUE),
      n_drop_flipper = sum(xb$filter_action == "drop_ambiguous_mirror_eqtl_artifact", na.rm = TRUE),
      genes = paste(genes_b, collapse = ";")
    )

    xb[, block_index := NULL]
    part_i <- part_i + 1L
    annotated_parts[[part_i]] <- xb
  }
}

annotated <- rbindlist(annotated_parts, use.names = TRUE, fill = TRUE)
annotated[, keep_after_filter := !filter_action %in% c(
  "drop_mhc",
  "drop_redundant_same_sign",
  "drop_ambiguous_mirror_eqtl_artifact",
  "drop_confirmed_mixed_artifact",
  "drop_possible_mixed_artifact"
)]
setcolorder(annotated, c(
  "gene", "disease", "fitness", "combo", "b1", "b1_q_fdr", "b2", "b2_q_fdr",
  "chr", "pos", "cm", "recomb_block_id", "recomb_block_n_genes",
  "recomb_block_disease_sign_pattern", "filter_action", "filter_reason",
  "representative_gene", "keep_after_filter"
))

after <- annotated[keep_after_filter == TRUE]
block_summary <- rbindlist(block_summaries, use.names = TRUE, fill = TRUE)
pair_evidence <- rbindlist(pair_evidence_rows, use.names = TRUE, fill = TRUE)

count_before <- combo[, .(before_n = .N), by = combo]
count_after <- after[, .(after_n = .N), by = combo]
count_summary <- merge(count_before, count_after, by = "combo", all = TRUE)
count_summary[is.na(before_n), before_n := 0L]
count_summary[is.na(after_n), after_n := 0L]
count_summary[, removed_n := before_n - after_n]
count_summary <- rbind(
  data.table(
    combo = "total",
    before_n = nrow(combo),
    after_n = nrow(after),
    removed_n = nrow(combo) - nrow(after)
  ),
  count_summary[order(combo)]
)
count_summary[, removed_fraction := removed_n / before_n]

count_by_fitness <- merge(
  combo[, .(before_n = .N), by = .(fitness, combo)],
  after[, .(after_n = .N), by = .(fitness, combo)],
  by = c("fitness", "combo"),
  all = TRUE
)
count_by_fitness[is.na(before_n), before_n := 0L]
count_by_fitness[is.na(after_n), after_n := 0L]
count_by_fitness[, removed_n := before_n - after_n]
count_by_fitness[, removed_fraction := removed_n / before_n]
count_by_fitness <- rbind(
  merge(
    combo[, .(before_n = .N), by = fitness],
    after[, .(after_n = .N), by = fitness],
    by = "fitness",
    all = TRUE
  )[, `:=`(combo = "total")][, .(fitness, combo, before_n, after_n)],
  count_by_fitness,
  fill = TRUE
)
count_by_fitness[is.na(before_n), before_n := 0L]
count_by_fitness[is.na(after_n), after_n := 0L]
count_by_fitness[, removed_n := before_n - after_n]
count_by_fitness[, removed_fraction := removed_n / before_n]

action_summary <- annotated[, .(
  combo_rows = .N,
  unique_genes = uniqueN(gene),
  kept_rows = sum(keep_after_filter),
  dropped_rows = sum(!keep_after_filter)
), by = filter_action][order(filter_action)]

action_by_fitness <- annotated[, .(
  combo_rows = .N,
  unique_genes = uniqueN(gene),
  kept_rows = sum(keep_after_filter),
  dropped_rows = sum(!keep_after_filter)
), by = .(fitness, filter_action)][order(fitness, filter_action)]

gene_action <- annotated[, .(
  combo_rows = .N,
  dropped_rows = sum(!keep_after_filter),
  actions = paste(sort(unique(filter_action)), collapse = ";"),
  any_drop = any(!keep_after_filter),
  all_drop = all(!keep_after_filter)
), by = gene][order(-any_drop, -all_drop, gene)]

manifest <- data.table(
  parameter = c(
    "significance_method",
    "significance_cutoff",
    "mixed_artifact_mode",
    "cm_gap_max",
    "mhc_exclusion",
    "shared_n_min",
    "shared_fraction_min",
    "same_sign_eqtl_cor_min",
    "flipper_eqtl_cor_max",
    "n_combo_rows_before",
    "n_combo_rows_after",
    "n_unique_genes_before",
    "n_unique_genes_after"
  ),
  value = c(
    significance_method,
    sig_cutoff,
    mixed_artifact_mode,
    cm_gap_max,
    "chr6:25000000-34000000",
    shared_n_min,
    shared_fraction_min,
    same_sign_eqtl_cor_min,
    flipper_eqtl_cor_max,
    nrow(combo),
    nrow(after),
    uniqueN(combo$gene),
    uniqueN(after$gene)
  )
)

fwrite(annotated[order(disease, fitness, chr, cm, gene)], file.path(out_dir, sprintf("%s_recomb_block_shared_iv_filter_combos_with_annotation.csv.gz", sig_tag)))
fwrite(after[order(disease, fitness, chr, cm, gene)], file.path(out_dir, sprintf("%s_recomb_block_shared_iv_filter_combos_after_filter.csv.gz", sig_tag)))
fwrite(block_summary[order(disease, fitness, chr, min_cm)], file.path(out_dir, sprintf("%s_recomb_block_summary.csv.gz", sig_tag)))
fwrite(pair_evidence[order(disease, fitness, recomb_block_id, gene_a, gene_b)], file.path(out_dir, sprintf("%s_recomb_block_pair_evidence.csv.gz", sig_tag)))
fwrite(count_summary, file.path(out_dir, sprintf("%s_combo_count_before_after_recomb_shared_iv_filter.csv", sig_tag)))
fwrite(count_by_fitness[order(fitness, combo)], file.path(out_dir, sprintf("%s_combo_count_before_after_recomb_shared_iv_filter_by_fitness.csv", sig_tag)))
fwrite(action_summary, file.path(out_dir, sprintf("%s_filter_action_summary.csv", sig_tag)))
fwrite(action_by_fitness, file.path(out_dir, sprintf("%s_filter_action_by_fitness.csv", sig_tag)))
fwrite(gene_action, file.path(out_dir, sprintf("%s_gene_action_summary.csv", sig_tag)))
fwrite(manifest, file.path(out_dir, "filter_manifest.csv"))

message("Done. Outputs written to: ", out_dir)
print(count_summary)
print(action_summary)
