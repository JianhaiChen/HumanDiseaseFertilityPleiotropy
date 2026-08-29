# SNP-level antagonism vs prevalence — bypassing FusioMR/FusioS entirely

Run 2026-08-19/20. Everything below uses **only the signs** of GWAS effect estimates.
No eQTL, no MR, no Bayesian estimate, no effect magnitude enters the analysis.

## Design

For each FinnGen R13 endpoint:
1. take SNPs with disease-side `P < threshold` (threshold only on the **disease** side;
   the fertility side is read unfiltered, which is what removes winner's curse);
2. LD-clump against the 1000G EUR panel (`r² = 0.1`, 1 Mb) to independent risk loci;
3. drop the extended MHC (6:25–34 Mb), 8p23.1 (8:8–12 Mb) and 17q21.31 (17:43–46 Mb)
   wholesale — r²-clumping does not clean long-range LD;
4. drop strand-ambiguous (A/T, C/G) and multi-rsID variants;
5. polarise each locus to its **risk allele** (β_D > 0);
6. read `sign(β_F)` for that same allele from the UK Biobank fertility GWAS
   (number of children fathered / mothered);
7. count **antagonistic** (risk allele raises fertility) vs **synergistic** (risk allele
   lowers fertility) loci.

Join key is the **rsID**: the fertility GWAS are hg19 and FinnGen R13 is hg38, so a
chrom:pos join is invalid (it loses ~98% of loci).

FinnGen × UK Biobank have no sample overlap, so β_D and β_F have independent errors.
The same analysis on the 97 UK Biobank diseases would not have this property.

## Headline result

Working point: disease-side `P < 1e-5`, endpoints with ≥50 independent loci
(**n = 485 endpoints, 119,833 loci**).

| | Pearson r vs log10(prevalence) | P |
|---|---|---|
| male fertility | **+0.271** | 1.3e-09 |
| female fertility | **+0.238** | 1.1e-07 |

Antagonistic share by prevalence quintile (male fertility):

| quintile | median prevalence | antagonistic | synergistic | share | binomial P |
|---|---|---|---|---|---|
| Q1 | 0.32% | 4,426 | 4,659 | **48.72%** | 1.5e-02 |
| Q2 | 0.98% | 6,471 | 6,583 | 49.57% | 0.33 |
| Q3 | 2.01% | 8,242 | 7,622 | 51.95% | 8.9e-07 |
| Q4 | 4.19% | 12,131 | 11,370 | 51.62% | 7.1e-07 |
| Q5 | 10.76% | 30,384 | 27,945 | **52.09%** | 5.7e-24 |

Antagonistic:synergistic ratio 0.950 (Q1) → 1.087 (Q5), Fisher P = 2.2e-09.
The sign crossover sits at roughly 1–2% prevalence: rarer diseases are
synergistic-dominated, commoner diseases antagonistic-dominated.

Interpretable effect size: **prevalence up one order of magnitude ⇒ antagonistic
share up ≈2 percentage points**.

## Controls

**Sign-permutation null** (flip the risk-allele assignment at random per locus,
everything else identical): share 49.99% (male) / 50.02% (female); correlation with
prevalence r = −0.031 / −0.042, both null.

**Positive control** — the SNP-level share tracks the existing model-based measures
(n_loci ≥ 100, non-WIDE, n = 252):

| | r | P |
|---|---|---|
| vs ρ_GE (male) | +0.499 | 3.0e-17 |
| vs r_g (male) | **+0.589** | 1.0e-24 |
| vs ρ_GE (female) | +0.480 | 6.0e-16 |
| vs r_g (female) | +0.523 | 5.5e-19 |

**Robustness grid** (all 9 combinations of n_loci ≥ 50/100/200 × all/non-WIDE/sex-shared):
male r = 0.25–0.36 (all P < 1e-4), female r = 0.21–0.27 (all P < 6e-3);
binomial GLM with ICD-chapter fixed effects, β = 0.066–0.117, all P < 2.3e-4.

**Covariate stress test** (`share ~ log10 prev + ICD chapter FE + onset age + mortality HR`):
male P = 6.9e-04, female P = 4.6e-03. Coefficients barely move.

## What the P-value ladder shows

Holding the **endpoint set fixed** (157 endpoints with ≥50 loci at 5e-8), tightening the
disease-side threshold monotonically *weakens* the signal:

| threshold | loci | share (male) | share (female) |
|---|---|---|---|
| 1e-4 | 201,693 | 51.17% | 51.06% |
| 1e-5 | 86,780 | 51.14% | 50.64% |
| 1e-6 | 48,382 | 50.81% | 50.39% |
| 1e-7 | 31,508 | 50.48% | 50.03% |
| 5e-8 | 28,324 | 50.35% | 49.90% |

The prevalence spread is essentially constant across these subsets (log10 SD 0.60–0.75),
so this is not a range-restriction artefact. **The antagonism lives in weak loci, not in
the strongest disease loci** — consistent with balancing selection maintaining
moderate-effect variation while large-effect risk alleles are either purged or purely
deleterious. Working range is therefore 1e-4 to 1e-6.

## Negative results

- **Ancestral/derived polarisation adds nothing.** Crossing the risk allele's derived
  state with the fertility sign: risk-derived 51.92% antagonistic vs risk-ancestral
  52.15% (male). No four-way structure beyond the two-way split.
- **Onset age**: gene-level antagonistic share moves +0.011/yr for male fertility
  (P = 0.011) and −0.016/yr for female fertility (P = 1.5e-04), i.e. under one
  percentage point across the whole 20→70 year range. The classic selection-shadow
  prediction (late onset ⇒ more antagonistic) is **not** supported; the female
  direction is opposite to it, and replicates at the SNP level (P = 1.7e-03).
- **Mortality (Risteys Cox HR)**: all associations with ρ_GE, r_g and antagonistic
  share vanish under ICD-chapter fixed effects (P = 0.63 for Δρ_GE). Pure category
  confounding.
- **Sex ratio of prevalence**: null throughout.

By contrast prevalence itself is completely unmoved by chapter fixed effects
(ρ_GE: P = 2.4e-81 → 1.5e-66 → 6.6e-43 with age and mortality added), so the
prevalence relation is not a proxy for disease category, severity, or onset timing.

## Files

- `snp_antag_prevalence.pdf/png` — panels: (a) scatter, (b) quintile deviation,
  (c) threshold ladder with permutation null
- `snp_antag_per_endpoint.tsv` — per-endpoint counts at all five thresholds
- `quintile_summary.tsv`, `threshold_ladder.tsv`

Cluster code: `randi:${BASEDIR}/selection/{snpantag.R,prep_fert.R,snpantag.sbatch}`,
outputs in `snpantag/chunks/`. Local analysis scripts in the session scratchpad.
