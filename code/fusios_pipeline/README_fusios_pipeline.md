# FusioS (FusioMR$_s$) ρ_GE pipeline — 给 A 整理 GitHub 用

FinnGen/UKB disease-fertility 表达效应相关 ρ_GE 的完整分析流程 + 代码文件位置。
所有代码在 **Randi** (`jianhaichen@randi.cri.uchicago.edu`)，scp 下来即可整理。

## 1. 方法概述
FusioMR$_s$ = single-exposure single-outcome (SESO) Gibbs：用 GTEx cis-eQTL 作 IV，
估计每个基因 genetically-regulated expression 对**单个** trait 的因果效应 (b_fusio, se_fusio)。
per-gene z = b_fusio/se_fusio；**ρ_GE = 跨基因 Pearson cor(disease_z, fertility_z)**。
disease 和 fertility 各自独立跑 FusioS（single-trait），再算相关。

## 2. 数据输入
| 数据 | 路径 (Randi) |
|---|---|
| FinnGen R13 disease GWAS (2217) | `/scratch/jianhaichen/finngen_r13/finngen_R13_<phenocode>.gz` |
| UKB disease GWAS (sex-specific std) | `/scratch/jianhaichen/fusios_muscle/sexspec_std/*.summaryfile` |
| father/mother fertility (UKB, 子女数) | `/gpfs/data/linchen-lab/Bowei/gwas_fitness/noffspring{father,mother}.fitness.summaryfile` |
| GTEx v10 Muscle cis-eQTL (rsid) | `/gpfs/data/linchen-lab/Bowei/dgtex/gtex_v10_rsid/eqtl_Muscle_Skeletal.v10.allpairs.chr{1-22}.txt.gz` |
| LD ref (per-chr bfile) | `/scratch/jianhaichen/eur_perchr/chr{1-22}.{bed,bim,fam}` |

Muscle_Skeletal 选它因为 GTEx 样本量最大。

## 3. IV 标准 (= main_july.tex Methods L146，全 pipeline 统一)
- cis-eQTL **P ≤ 0.001**
- LD-clump：PLINK，**100-kb window，r² = 0.1**
- 每基因 **≥ 5** 个独立 IV（否则排除）
- top-N by eQTL p（FinnGen 30 / UKB 50；纯计算优化，near-lossless）

## 4. Pipeline（每 gene × chr）
1. 读 eQTL (p<0.001) + trait GWAS → 按 rsid harmonize，align 到 eQTL effect allele
2. ≥5 IV floor + top-N by p
3. per-gene clump（**per-chr bfile**）→ SESO Gibbs (`process_gene_bk_only`)
4. 输出 per-gene：`gene, niv, b_fusio, se_fusio, p_fusio`

## 5. 代码文件 (Randi — 要放 GitHub 的)
| 文件 | 作用 |
|---|---|
| `/scratch/jianhaichen/fusios_finngen/code/fusios_finngen.R` | FinnGen disease FusioS 主脚本 (22 chr in one job) |
| `/scratch/jianhaichen/fusios_muscle/code/fusios_ukb.R` | UKB disease FusioS 主脚本 |
| `/scratch/jianhaichen/fusios_finngen/code/functions_perchr.R` | clump 函数副本，bfile→全局 `PERCHR_BFILE`（per-chr 加速用） |
| `/scratch/jianhaichen/fusios_finngen/agg_rhoGE.R` | 聚合 per-gene z + 算 ρ_GE (vs father/mother) |
| `/scratch/jianhaichen/fusios_finngen/fusios_finngen.sbatch` | SLURM array 提交脚本 |
| **FusioMR core (Bowei，依赖)** `/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/code/` | `functions.R`, `init_setup_seso.R`, `gibbs_seso_uhp_only.cpp` (SESO Gibbs C++) |

`functions_perchr.R` = 把 Bowei 的 `functions.R` 里 6 处硬编码 bfile
(`education_AD_GWAS/EUR`) 替换成全局变量 `PERCHR_BFILE`；主脚本每条 chr 循环里设
`PERCHR_BFILE <- eur_perchr/chr<N>`。per-chr bfile 由
`plink --bfile EUR --chr N --make-bed` 生成，byte-identical，**~7× 加速**。

## 6. ρ_GE 计算 (`agg_rhoGE.R`)
- per gene z = b_fusio/se_fusio，**|z| < 9** cap
- ρ_GE = Pearson cor(disease_z, father_z) 和 (disease_z, mother_z) over shared genes
- **≥ 30** shared genes；**MHC (chr6:25–34 Mb) 排除**（canonical 版需加，见注）
- 输出 `rhoGE_full.tsv`：`phenocode, rho_father, ng_f, rho_mother, ng_m`

## 7. 依赖
R 包：`data.table, tidyverse, pbapply, Rcpp`；PLINK 1.9；GTEx v10 rsid-annotated eQTL。
环境：`module load gcc; module load R/4.2.1`。

## 8. 计算注意（写 methods/README 可提）
瓶颈是**文件系统 I/O**（每基因 clump 读 bfile + eQTL），不是 CPU。
SLURM 并发 sweet spot ~186；>~470 并发会 I/O 饱和、吞吐崩溃。
per-chr bfile 是关键加速（全基因组 850万 variant → 每 chr ~85万）。
