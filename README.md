# Why are common diseases common?

My background is in evolutionary genetics, and my earlier work was on
evolutionarily new genes:

> Chen J-H, Landback P, Arsala D, Guzzetta A, Xia S, ... Zhang YE, Cheng J,
> Shen B, Long M. Evolutionarily new genes in humans with disease phenotypes
> reveal functional enrichment patterns shaped by adaptive innovation and
> sexual selection. *Genome Research* 35:379-392 (2025).
> doi:10.1101/gr.279498.124

That study used rare diseases, as did the variant-screening work I did
alongside it. One of those was a DMD splice-donor mutation causing Duchenne
muscular dystrophy with cryptorchidism:

> Chen J, Jia Y, Zhong J, Zhang K, Dai H, He G, Li F, Zeng L, Fan C, Xu H.
> Novel mutation leading to splice donor loss in a conserved site of DMD gene
> causes Duchenne muscular dystrophy with cryptorchidism.
> *Journal of Medical Genetics* 61:741-749 (2024). doi:10.1136/jmg-2024-109896

Duchenne affects about 1 in 5,000 boys. Every screening paper opens with a
number like that, and after enough of them I started wondering what prevalence
means in evolutionary terms. Duchenne is rare and stays rare: it is lethal
before reproduction, so selection removes it every generation. Type 2 diabetes
is not rare. Something has to be holding the alleles in.

I began working on this in the Lin lab. The first pass covered 97 diseases and
found nothing strong, because there was no prevalence data to work with. But
the intuition held: a common disease is common because its alleles keep being
transmitted, and transmission is fitness. That is an evolutionary question
before it is an epidemiological one.

So I added reproductive data. That is where prevalence turned out to track
antagonistic and synergistic pleiotropy between disease and fertility.

This repository holds the code for that analysis.

---

Code for the disease-fertility pleiotropy paper. Directories follow the Methods sections.

FusioMR_m = joint model (97 UK Biobank diseases, 50 GTEx tissues).
FusioMR_s = single-trait model (FinnGen layer; also the shared-instrument check in muscle).

00_pipeline_templates/
  preclump_R_full.R + preclump_full.sbatch    instrument construction
  run_disease_full.sbatch, run_dz_chunk.sbatch, run_fert_full.sbatch
                                              transcriptome-wide MR submission
  acat50_full.sbatch, acat_fert_full.sbatch   gene-level ACAT submission
  README.md                                   stage order and the settings that
                                              must not change mid-campaign

01_transcriptome_wide_mr/
  run_fusiom_fast.R           FusioMR_m runs
  gene_level_acat_summary.R   gene-level ACAT -> b1 (disease) / b2 (fertility) tables
  acat_female.R               compact ACAT reference; signed effect = min-P tissue

02_local_filtering/
  apply_recomb_block_shared_iv_filter_fdr005.R
    Paper used SIGNIFICANCE_METHOD=bonf_global SIG_CUTOFF=0.05 MIXED_ARTIFACT_MODE=broad
    (outputs named bonf050_broadmixed_*), despite the filename.

03_coupling/
  fig_BC_97_bysex.R     per-disease rho_GE, sex-specific
  validation_sex.R      the two rho_GE sensitivity analyses
  mhc_compare.R         MHC in/out
  ldsc/                 LDSC r_g pipeline (clean -> align -> munge -> rg)
  REPRODUCE_ldsc_and_rhoGE.md

04_prevalence_finngen/
  proxy_prev.R              FinRegistry prevalence vs FinnGen case fraction
  fusios_iv5.R              50-tissue FinnGen run
  acat50.R                  ACAT across tissues, sex-incompatible tissues dropped
  noise_deconvolution.R     false-positive correction; runs as-is off ../data/
                            (cut from the manuscript, kept for review response)

05_sex_limited/
  classify_endpoints.R          four-class rules; fertility-bound is ICD chapter XV only
  endpoint_sex.tsv
  delta_rg_test.R               delta r_g per endpoint
  make_figS_maletract.R         includes the Steiger test
  make_tableS9_sexconflict.R    Hodges-Lehmann + Wilcoxon interval

data/   small tables so the above can be checked without the cluster

Not included: deconv_antag.R, superseded by noise_deconvolution.R.

06_figures/
  panel and figure scripts
  FIGURES.md          each figure -> the script behind it
  assemble_fig5.py    panel placement step

07_fusiomr_core/
  FusioMR Gibbs samplers and setup, sourced by the two MR drivers

Paths are variables throughout: ${ROOT} ${PROJ} ${DATA} ${MR} ${BASEDIR}
${LABDIR} ${FUSIOMR_CODE} ${REFPANEL}. Set them before running.
