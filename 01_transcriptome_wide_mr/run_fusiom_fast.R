#!/usr/bin/env Rscript
## FusioMR_m transcriptome-wide MR. One job per (tissue, disease, fitness)
## triple, all 22 chromosomes inside the job: the local prior is read once, the
## global prior built once, the Rcpp core compiled once.
##
## Results accumulate in memory and are written to one file per job. The Rcpp
## cache is read from a shared prebuilt copy and rebuilt node-locally, since a
## shared cache directory is not safe across concurrent array tasks.
## Genes are subset by data.table key and run under mclapply.

suppressPackageStartupMessages({
  library(data.table); library(mvtnorm); library(LaplacesDemon)
  library(dirmult); library(invgamma); library(parallel)
})

source("${FUSIOMR_CODE}/conflict/code/set_variance_prior.R")
source("${FUSIOMR_CODE}/conflict/code/init_setup_semo.R")

## ---- hyperparameters: identical to run_fusiom.R lines 19-20 ----
glob_q <- 0.75; hybrid <- TRUE; kappa_hybrid <- 50
c_gamma <- 4; c_theta <- c_gamma * 1.5
kappa_gamma_m <- 1; kappa_theta_m <- 2
niter <- 20000

## Reproducibility: the Gibbs sampler has no internal seed. L'Ecuyer-CMRG is
## required for reproducible streams across mclapply forks (set.seed alone is not).
RNGkind("L'Ecuyer-CMRG"); SEED <- 20260807L

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: run_fusiom_fast.R <tissue> <disease> <fitness>")
tissue <- args[1]; disease <- args[2]; fitness <- args[3]

BASE    <- "${BASEDIR}/rerun_conflict"
OUTDIR  <- Sys.getenv("FUSIOM_OUTDIR", file.path(BASE, "mr_res_fast"))
NCORE   <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(OUTDIR, paste0(tissue, "__", disease, "__", fitness, ".mrres.txt"))
if (file.exists(outfile) && file.info(outfile)$size > 0) {
  cat("Already done, skipping:", outfile, "\n"); quit(save = "no")
}

## ---- Rcpp: shared prebuilt cache, node-local writes ----
## sourceCpp rewrites its cache index on every call, so pointing hundreds of
## concurrent array tasks at one shared cacheDir corrupts it (load() fails with
## "magic number"). Read the shared copy, write to node-local scratch.
CC_SRC <- "${BASEDIR}/rerun_conflict/ccache_shared"
CC <- file.path(Sys.getenv("TMPDIR", "/tmp"),
                paste0("fmcache_", Sys.getenv("SLURM_JOB_ID", "x"), "_",
                       Sys.getenv("SLURM_ARRAY_TASK_ID", "0")))
dir.create(CC, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(CC_SRC))
  file.copy(list.files(CC_SRC, full.names = TRUE), CC, recursive = TRUE, overwrite = TRUE)
Rcpp::sourceCpp("${FUSIOMR_CODE}/conflict/code/gibbs_semo_uhp_only_bound.cpp",
                cacheDir = CC)
on.exit(unlink(CC, recursive = TRUE), add = TRUE)

## ---- global priors: built once for the whole trait pair ----
build_global_mom <- function(local_mom_vec, q = 0.5, trim = 0.05, min_pos = 20) {
  x <- local_mom_vec[is.finite(local_mom_vec) & local_mom_vec > 0]
  if (length(x) == 0) return(0)
  if (length(x) < min_pos) return(stats::median(x, na.rm = TRUE))
  lo <- stats::quantile(x, trim, na.rm = TRUE); hi <- stats::quantile(x, 1 - trim, na.rm = TRUE)
  x <- x[x >= lo & x <= hi]
  if (!length(x)) return(0)
  stats::quantile(x, q, na.rm = TRUE)
}
build_global_mom_covmat <- function(mat, q = 0.5, eps = 1e-8) {
  A <- pmax(mat[, 1], eps); B <- pmax(mat[, 2], eps); C <- mat[, 3]
  r <- pmin(pmax(C / sqrt(A * B), -1 + 1e-10), 1 - 1e-10)
  s1 <- exp(stats::quantile(log(A), q, na.rm = TRUE))
  s2 <- exp(stats::quantile(log(B), q, na.rm = TRUE))
  r_bar <- tanh(stats::quantile(atanh(r), 0.5, na.rm = TRUE))
  matrix(c(s1, r_bar * sqrt(s1 * s2), r_bar * sqrt(s1 * s2), s2), 2, 2)
}

lp_file <- file.path(BASE, "local_priors", paste0(tissue, "__", disease, "__", fitness, ".local_prior.txt"))
if (!file.exists(lp_file) || file.info(lp_file)$size == 0) {
  cat("Missing/empty local prior file. Exit.\n"); quit(save = "no")
}
lp <- fread(lp_file)
setnames(lp, c("tissue","disease","fitness","chr","gene","K0","K1",
               "mom_gamma","mom_theta11","mom_theta22","mom_theta12","beta01","beta02"))
lp <- unique(lp, by = "gene")

global_m_gamma <- build_global_mom(lp$mom_gamma, q = glob_q, trim = 0.05, min_pos = 20)
global_m_theta_mat <- build_global_mom_covmat(as.matrix(lp[, mom_theta11:mom_theta12]),
                                              q = glob_q, eps = 1e-8)
setkey(lp, chr)

## ---- per-gene worker ----
run_gene <- function(df, gg, chr) {
  tryCatch({
    b_exp <- df$beta_eqtl; se_exp <- df$se_eqtl
    b1v <- df$beta_d_h; se1v <- df$se_d
    b2v <- df$beta_f_h; se2v <- df$se_f
    K <- nrow(df)

    vp <- set_variance_priors_m2(
      ghat = b_exp, gse = se_exp,
      Ghat_mat = matrix(c(b1v, b2v), ncol = 2), Gse_mat = matrix(c(se1v, se2v), ncol = 2),
      beta0 = NULL, K = K, Kmin = 5, Kmax = 20, rho12 = 0, rho1g = 0, rho2g = 0,
      c_gamma = c_gamma, c_theta = c_theta,
      global_mean_gamma = global_m_gamma, global_mean_theta = global_m_theta_mat,
      hybrid = hybrid, kappa_hybrid = kappa_hybrid, z_thresh = NULL, trim = 0.0,
      kappa_gamma = kappa_gamma_m, kappa_theta = kappa_theta_m)

    sv <- init_setup_semo_uhp_only(niter, K, beta_1_init = vp$beta0[1],
                                   beta_2_init = vp$beta0[2],
                                   sigma_gamma_init = sqrt(vp$gamma$prior_mean))

    r <- gibbs_semo_uhp_only_rcpp(niter, K,
      sv$beta_1_tk, sv$beta_2_tk, sv$theta_1_tk, sv$theta_2_tk,
      sv$gamma_tk, sv$sigma2_gamma_tk, b1v, b2v, se1v^2, se2v^2, b_exp, se_exp^2,
      vp$gamma$a, vp$gamma$b, vp$theta$prior_mean, vp$theta$nu, vp$theta$Phi)

    ids <- (niter / 2 + 1):niter
    b1 <- mean(r$beta_1_tk[ids], na.rm = TRUE); s1 <- sd(r$beta_1_tk[ids], na.rm = TRUE)
    b2 <- mean(r$beta_2_tk[ids], na.rm = TRUE); s2 <- sd(r$beta_2_tk[ids], na.rm = TRUE)

    data.table(tissue = tissue, disease = disease, fitness = fitness, chr = chr, gene = gg,
      K = K,
      b1_semo = b1, se1_semo = s1,
      p1_semo = 2 * exp(pnorm(abs(b1) / s1, lower.tail = FALSE, log.p = TRUE)),
      b2_semo = b2, se2_semo = s2,
      p2_semo = 2 * exp(pnorm(abs(b2) / s2, lower.tail = FALSE, log.p = TRUE)))
  }, error = function(e) { cat("ERROR gene", gg, ":", conditionMessage(e), "\n"); NULL })
}

## ---- loop chromosomes, parallelise genes ----
ACC <- list()
## loop variable is cc, not chr: `chr` is a column of lp, and using it as the
## key-lookup value would be ambiguous inside lp[...]
for (cc in 1:22) {
  f <- file.path(BASE, "df_clumped",
                 paste0(tissue, "__", disease, "__", fitness, "__chr", cc, ".txt.gz"))
  if (!file.exists(f) || file.info(f)$size == 0) next
  dfc <- fread(f)
  need <- c("gene","beta_eqtl","se_eqtl","beta_d_h","se_d","beta_f_h","se_f")
  if (length(setdiff(need, names(dfc)))) stop("Missing columns in ", f)

  genes <- intersect(unique(dfc$gene), lp[.(cc), gene, nomatch = NULL])
  if (!length(genes)) next
  cat(tissue, disease, fitness, "chr", cc, ":", length(genes), "genes\n")

  ## split() once beats filtering the whole table per gene
  parts <- split(dfc[gene %in% genes], by = "gene", keep.by = TRUE)

  set.seed(SEED + cc)
  res <- mclapply(genes, function(g) run_gene(parts[[g]], g, cc),
                  mc.cores = NCORE, mc.preschedule = TRUE)
  ## run_gene returns a data.table or NULL; anything else is a worker crash
  ok <- vapply(res, is.data.table, logical(1))
  if (sum(!ok)) cat("  ", sum(!ok), "genes failed on chr", cc, "\n")
  if (any(ok)) ACC[[length(ACC) + 1L]] <- rbindlist(res[ok])
}

final <- if (length(ACC)) rbindlist(ACC, fill = TRUE) else data.table()
fwrite(final, outfile, sep = "\t")
cat("Wrote", nrow(final), "rows ->", outfile, "\n")
