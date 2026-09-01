# Repository contents

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

--------------------------------------------------------------------------
Added 2026-09-01, for the Nature Communications version
--------------------------------------------------------------------------

04_prevalence_finngen/joint_model/
  joint_model_prevalence.R    the single pre-specified model behind the
                              "beyond mortality and disease onset" section.
                              Every partial R^2 in that section and in
                              Fig. 4d-f comes from this one fit; the model is
                              never re-specified per term. Prints the numbers
                              quoted in the text (coupling 0.206, mortality
                              0.152, onset 0.047, spline-over-linear
                              dR^2 = 0.028) and the onset-tertile description
                              used in Fig. 4f.
  make_supp_table12.R         builds Supplementary Table 12 from that fit.
  model_input_1678.tsv        its input: coupling, log 5-year mortality HR,
                              median age at first recorded event and ICD
                              chapter for the 1,678 endpoints with a Risteys
                              mortality estimate.

  NOTE on excluded covariates: case count, locus count and the number of genes
  behind each coupling estimate are all downstream of the disease case count,
  which is nearly collinear with prevalence (r = 0.97, log10 scale). Adjusting
  for any of them removes the outcome variance the model is asked to explain.

04_prevalence_finngen/
  rebuild_table10_fdr05.R     rebuilds Supplementary Table 10 (directional
                              classes) at FDR < 0.05.

06_figures/rebuilt_2026-09/
  rebuild_class_figures_fdr05.R   directional-class figures at Benjamini-
                                  Hochberg FDR < 0.05 (the manuscript switched
                                  from nominal P < 0.05 on 2026-08-31).
                                  Rebuilds main Fig. 4a and Supplementary
                                  Figs. 10 and 13.
  fig4_panels_de_partial_residual.R  Fig. 4d (partial-residual plot for the
                                  onset spline) and Fig. 4e (added-variable
                                  plots for the two coupling terms).
  replot_fig5b_recomb.R           Fig. 5b redrawn from the cached
                                  recomb_by_group.tsv; the original script
                                  reads "Supplementary Table_v6.xlsx", which
                                  has since been renamed.
  fig4_6panel.tex                 six-panel Fig. 4 assembly
  fig5_4panel_wide.tex            four-panel Fig. 5 assembly
  fig1_panelA_hypotheses.tex      alternative Fig. 1a (evolutionary model
                                  schematic); not used in the submitted
                                  version, which keeps the analysis-framework
                                  panel.
