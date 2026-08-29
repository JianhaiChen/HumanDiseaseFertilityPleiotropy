# Reproducibility — disease–fertility genetic correlation (97 diseases) & Fig 3

Last updated 2026-07-12. Everything below reproduces the finalized main-text **Fig 3
(`Fig_BC_final.pdf`)** and the 97-disease disease–fertility correlation tables.
All scripts are in `./scripts/`. Two machines: **local Mac** (figures, ρ_GE) and
**Randi HPC** (LDSC r_g); `ssh randi`.

---

## 0. Key outputs (this directory)
| File | What |
|---|---|
| `dz_rg_97_ALLESTIMABLE.txt` | **FINAL** LDSC r_g, 97 diseases, 0 NA, sex-specific father(M)/mother(F)/lifespan |
| `core_rg_table.csv` | per-disease r_g + ρ_GE, sign/sig/#genes |
| `why_ns_table.csv` | non-significance diagnosis (n_strong, sex-discord) |
| `mhc_compare_table.csv` | ρ_GE with vs without MHC |
| `rhoge_auth_compare.csv` | FusioMR vs FusioS vs LDSC ρ_GE |
| `../Fig_BC_final.pdf`, `fusio2trait/Fig_BC_final.pdf` | **main-text Fig 3** (97 dz, sex-split) |

`dz_rg_97_ALLESTIMABLE.txt` **supersedes** the older `dz_rg_all.txt` (87 dz) and
`dz_rg_97_final.txt` (had 3 NA).

---

## 1. LDSC r_g — 97 diseases, sex-specific (on Randi)
Working dir: `${BASEDIR}/selection/`. Env: `conda activate ldsc` (python 2.7).
Reference: `eur_w_ld_chr/` (EUR LD scores) + `w_hm3.snplist`.
Fitness sumstats already munged: `father.al.sumstats.gz`, `mother.al.sumstats.gz`, `lifespan.al.sumstats.gz`.

Pipeline per disease: **raw GWAS → clean → align to HM3 → munge → rg vs father/mother/lifespan**
- `process_one.sh dz file N` = build_clean2.sh → align_one.sh → munge_sumstats.py (original 87)
- `clean6_fix.awk` = **fixed cleaner** (handles the formats that broke build_clean2.sh):
  - GWAS-Catalog harmonised with NA `beta` → fall back to `log(odds_ratio)`
  - PGC format: SNP col = `ID`, p col = `PVAL`
  - FinnGen format: lowercase `alt`/`ref` alleles, SNP col = `rsids`
- munge uses `ldsc/munge_sumstats.py` (NOT ldsc.py); rg uses `ldsc/ldsc.py --rg`

### 1a. The 87 baseline + 7 recovered
- Original 87 valid from prior run (`dz_rg_all.txt`).
- `run6c.sh` — recovered 6 (Cataract, Melanoma, Pancreatic_Cancer, Panic, non-Hodgkins_lymphoma, thyroidcancer): raw GWAS existed in `full_manifest_final.tsv` but failed cleaning; clean6_fix.awk fixes them.
- autism — h2=0.20 was fine; just re-ran rg (father=-0.21).

### 1b. The 3 FinnGen rescues (Randi has internet; wget)
Original GWAS for these had **Mean Chi2 < 1 = no detectable SNP-h2** (LDSC out of bounds).
Replaced with better-powered **FinnGen R12** (European registry, full sumstats, all Mean Chi2>1):
- `finngen_brain.sh`: braincancer ← `finngen_R12_C3_BRAIN_EXALLC` (Malignant neoplasm of brain, 1894 cases)
- `fg_two.sh`: Sepsis ← `AB1_OTHER_SEPSIS` (17133); Vitamin_D ← `E4_VIT_D_DEF` (614)
- URL pattern: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_<PHENO>.gz`
- N passed = Ncase+Nctrl; sanity check: Sepsis×lifespan rg=-0.74 (sepsis shortens life ✓)

### 1c. Assemble final file
On Randi, replace the 3 NA rows in `dz_rg_97_final.txt` with the FinnGen rg values →
`dz_rg_97_ALLESTIMABLE.txt` (97 rows, 0 NA). scp to this dir + mr_res root + session dir.

**h2-out-of-bounds cause** (diagnosed, NOT outliers): Mean Chi2<1 / Lambda GC<1 =
deflated test statistics = GWAS has no polygenic signal (underpowered / rare disease).
Fixed only by a better-powered GWAS, not by method tweaks.

---

## 2. ρ_GE (gene-expression correlation) — local, all 97
Authoritative method (in `fig_BC_97_bysex.R`, matches original `fig_BC_final.R`):
```
z1 = disease top_beta/top_se  (gene_level_b1_disease_acat_summary.dropbox.csv)
z2 = fitness top_beta/top_se  (gene_level_b2_fitness_acat_summary.dropbox.csv)
exclude MHC (chr6:25-34Mb); cap |z|<9; per disease×sex: rge = cor(z1,z2), ng>=30
```
ρ_GE is per-sex (male/female) — **never averaged** (33 diseases are sex-antagonistic;
averaging cancels real signal). See `mhc_compare.R`/`mhc_only.R` for the MHC decision:
genome-wide ρ_GE robust to MHC (r=0.94); within-MHC is one sex-flipping haplotype → keep MHC out.

---

## 3. Main-text Fig 3 — `scripts/fig_BC_97_bysex.R`
```
cd ${MR}/session_2026-07-09_iv_filter
Rscript fig_BC_97_bysex.R
```
Inputs (all STABLE, no ephemeral scratchpad): the two gene_level ACAT csvs,
`gene_coords_v19.tsv`, `dz_rg_97_ALLESTIMABLE.txt`, `Supplementary_Table_1.xlsx`,
`Supplementary_Table_1_corrected_v2.csv` (in mr_res/), `recomb/` maps (mr_res/recomb/),
and `colleague_shared_iv_filter.csv` (session_2026-07-09_iv_filter/, also copied here).
Paths are set at the top of the script (SC=mr_res for recomb; CFPATH=colleague file) — self-contained.
Outputs: `fusio2trait/Fig_BC_final.pdf` (main text) + `Fig3_BC_97_bysex.pdf` (alias) + mr_res root copy.

- **Panel A**: disease-fertility correlation vs % antagonistic loci, faceted by method
  (ρ_GE MR / r_g LDSC), colored by sex. All 97 (threshold n>=1); points sized by #independent
  loci; regression AND reported r/P are **n-weighted** (cov.wt r + weighted-lm p) so low-n
  diseases neither dilute nor inflate. → ρ_GE M r=0.51/F r=0.55; r_g M r=0.26/F r=0.38.
- **Panel B**: per-disease r_g & ρ_GE lollipops, **faceted by sex (Male/Female), not averaged**.

---

## 4. One-command re-run order
1. (Randi) rerun 1a–1c if raw GWAS change → new `dz_rg_97_ALLESTIMABLE.txt`, scp to mr_res root.
2. (local) `Rscript fig_BC_97_bysex.R` → Fig_BC_final.pdf.
All inputs are on stable disks (no ephemeral temp dirs). `colleague_shared_iv_filter.csv`
and all analysis scripts are also copied into this dir (`./` and `./scripts/`) for self-containment.
