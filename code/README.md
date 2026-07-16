# Figure 3 — disease–fertility genetic correlation (code)

`fig3_genetic_correlation.R` computes the per-disease disease–fertility genetic
correlation with two estimators, **sex-specifically** (father=male, mother=female),
and writes the Figure-3 source table `fig3_disease_correlations.csv`.

| estimator | what | source |
|---|---|---|
| `rhoGE_male/female` | gene-expression correlation `cor(z_disease, z_fitness)` across genes, z = FusioMR MR effect (top_beta/top_se) | `gene_level_b1/b2_..._acat_summary` |
| `rgM` (father), `rgF` (mother) + SE | genome-wide SNP genetic correlation, LD-score regression | `dz_rg_97_ALLESTIMABLE.txt` |

## MHC handling — the one thing to get right
The two estimators must treat the MHC region (chr6:25–34 Mb) **identically**, and
here **both EXCLUDE it** (`INCLUDE_MHC <- FALSE`).

- **r_g** uses the standard European LD-score reference `eur_w_ld_chr`, which has
  **no SNPs in chr6:24,999,740–34,005,069** (a 9-Mb gap = the MHC). `ldsc.py --rg`
  can only use SNPs present in both the sumstats and the reference, so MHC SNPs are
  dropped automatically → **our r_g already EXCLUDES MHC** (verified directly from
  `eur_w_ld_chr/6.l2.ldscore.gz`).
- `rho_GE` therefore also excludes MHC, via the explicit `!isMHC(chr,pos)` filter,
  so both measures are on the same footing.
- **Rationale:** the MHC's extreme long-range LD prevents reliable causal-gene
  identification in gene-level MR/TWAS and biases SNP-level LDSC; it is removed as
  a region. MHC genes are still retained (flagged `in_MHC`) in Supplementary
  Table 2 — exclusion is a causal-resolution limitation, not a biological claim.
- **Robustness:** for `rho_GE`, including MHC vs treating it as ~7 recombination-
  block representatives changes per-disease `rho_GE` by ≤0.02 (male cor 1.00,
  female 1.00) — negligible.

## Result (both estimators exclude MHC — consistent)
- rho_GE vs r_g convergence: **male (father) r=0.68 P=2.5e-14, 72% sign-concordant;
  female (mother) r=0.70 P=1.0e-15, 82% sign-concordant.**
- 97 diseases, both estimators.

## Run
```
Rscript fig3_genetic_correlation.R    # edit DATA path at top first
```
Writes `fig3_disease_correlations.csv` (Fig 3 source data) and prints convergence stats.
