# Conflict effects of gene expression on complex diseases and fitness

Reproducibility code and key data tables for the disease–fertility pleiotropy /
antagonistic-pleiotropy analyses. Chen, Kang, Long & Chen.

The pipeline uses **FusioMR** (transcriptome-wide Mendelian randomization,
GTEx v10 eQTL × disease/fertility GWAS) to estimate, per gene, the effect of
genetically regulated expression on 97 diseases and on male (father) and female
(mother) fertility, and characterizes disease–fertility pleiotropy and its
relationship to a sex-specific disease–fertility genetic correlation.

## Repository layout
```
code/    analysis + figure scripts (R)
data/    small, shareable per-figure source data
```
Figures are produced by the scripts in `code/` (not committed); the manuscript
supplementary tables are distributed with the paper, not in this repo.

## Data provenance (large source data — deposited on Zenodo)
These files are too large for GitHub and will be deposited on **Zenodo
(DOI: TO_BE_ADDED)**. Download them and point the `DATA` variable at the top of
each script to the download folder.
| data | what | size |
|---|---|---|
| `gene_level_b1_disease_acat_summary.dropbox.csv` | per-gene disease MR (ACAT across tissues) | 937 MB |
| `gene_level_b2_fitness_acat_summary.dropbox.csv` | per-gene fertility MR (father/mother) | 20 MB |
| `fusios_independent_iv/agg_traits/*.txt` | FusioS per-gene single-trait MR (Muscle_Skeletal), 99 traits | 295 MB |
| `fusios_intermediate.tar.gz` | FusioS raw per-chromosome outputs + run scripts | 90 MB |
| `gene_coords_v19.tsv` | gene → chr,pos (hg19) | small |

LDSC munge+rg pipeline (Randi `/scratch/jianhaichen/selection/`) produced
`dz_rg_97_ALLESTIMABLE.txt`, which IS included in `data/`.
Small derived tables needed to reproduce the figures are in `data/`.

## MHC handling (consistent across ALL analyses)
The MHC region (chr6:25–34 Mb) is **excluded from every genetic-correlation and
pleiotropy analysis**, because its extreme long-range LD prevents reliable
causal-gene identification in gene-level MR/TWAS and biases SNP-level LDSC.
- r_g: the LDSC reference `eur_w_ld_chr` has a 9-Mb SNP gap over the MHC, so
  LDSC drops it automatically (verified).
- ρ_GE and all gene-level analyses: explicit `!isMHC(chr,pos)` filter.
- MHC genes are retained (flagged `in_MHC`) in Supplementary Table 2 (with the paper).

## Figure → code → data map
| figure | code | source data |
|---|---|---|
| **Fig 3** disease–fertility correlation & antagonism | `code/fig3_genetic_correlation.R` (data) + `code/fig3_plot_panels.R` (panels A/B/C) | `code/fig3_disease_correlations.csv`, `data/dz_rg_97_ALLESTIMABLE.txt`, gene_level b1/b2 |
| **SuppFig 6** % genes per 4-class vs correlation | `code/suppfig6_4class.R` | gene_level b1/b2, colleague IV-filter |
| **SuppFig 7** SNP r_g predicts antagonism / **SuppFig 8** ρ_GE robust to shared IV | `code/suppfig7_8_validation.R` | `data/rhoge_auth_compare.csv`, gene_level b1/b2 |
| **SuppFig** male vs female disease–fertility correlation (r_g \| ρ_GE) | `code/suppfig_sex_concordance.R` | `code/fig3_disease_correlations.csv`, Supplementary Table 1 (categories) |
| FusioMR vs FusioS ρ_GE (authoritative, MHC-excluded) | `code/fusioMR_vs_fusioS_rhoGE.R` | FusioS agg + gene_level b1/b2 |
| (exploratory, not in manuscript) sexual antagonism | `code/exploratory_sexual_antagonism.R` | FusioS agg father/mother |

### Not yet packaged (authoritative scripts to be added)
Fig 1 (study design / pleiotropy overview), Fig 2 (4-class antagonism),
Fig 4 (parent→daughter evolution), Fig 5 (tissue intralocus switching), and
SuppFig 1–5 (disease-pleiotropy / 4-class distributions). Their data come from
the same gene_level MR tables and Supplementary Tables 2–4.

## Key data tables (`data/`)
- `dz_rg_97_ALLESTIMABLE.txt` — LDSC r_g, 97 diseases, sex-specific (father/mother).
- `rhoge_auth_compare.csv` — FusioMR vs FusioS ρ_GE per disease×sex.
- `code/fig3_disease_correlations.csv` — Fig 3 source table (ρ_GE + r_g, 97 dz).
