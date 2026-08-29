# Why are common diseases common?

Human diseases vary enormously in prevalence. Some affect a substantial
fraction of the population, whereas others are exceedingly rare. Epidemiology
has documented this variation in extraordinary detail, yet a basic evolutionary
question remains:

**Why are some diseases common while others are rare?**

## Where this question came from

My background is in evolutionary genetics. Earlier in my career I worked in a
*Drosophila* laboratory, where I was trained in a tradition that places
enormous value on experimental precision: perturb a gene, control the genetic
background, observe the phenotype, and use reverse genetics to establish
mechanism.

That approach is extraordinarily powerful. But I became increasingly
uncomfortable with a view I sometimes encountered, that human genetic studies
are inherently less rigorous or less informative because humans cannot be
experimentally manipulated with the precision possible in model organisms.

My experience in human rare-disease genomics led me to a different view.

Human biology comes with a different kind of experimental resource: an enormous
natural record accumulated through centuries of medicine and, more recently,
through population-scale genomics. Human phenotypes have been observed across
diseases, organs, ages, environments, populations, and generations at a scale
and breadth that no experimental model can reproduce.

Consider Duchenne muscular dystrophy. Generations of physicians, geneticists,
pathologists, molecular biologists, and families have collectively
characterized its incidence, inheritance, natural history, reproductive
consequences, molecular pathology, and clinical progression. It would be
extraordinarily difficult to reproduce that depth of accumulated knowledge for
a single analogous phenotype in a model organism, let alone for thousands of
phenotypes.

This difference shaped how I came to think about human genetics. Experimental
control is one form of scientific power; accumulated variation and observation
at population scale is another.

My earlier work examined how evolutionarily new genes contribute to human
disease:

> Chen J-H, Landback P, Arsala D, Guzzetta A, Xia S, ... Zhang YE, Cheng J,
> Shen B, Long M. Evolutionarily new genes in humans with disease phenotypes
> reveal functional enrichment patterns shaped by adaptive innovation and
> sexual selection. *Genome Research* 35:379-392 (2025).
> doi:10.1101/gr.279498.124

In parallel, I worked directly on rare Mendelian disease, including the
identification of a splice-donor mutation in *DMD* causing Duchenne muscular
dystrophy:

> Chen J, Jia Y, Zhong J, Zhang K, Dai H, He G, Li F, Zeng L, Fan C, Xu H.
> Novel mutation leading to splice donor loss in a conserved site of DMD gene
> causes Duchenne muscular dystrophy with cryptorchidism.
> *Journal of Medical Genetics* 61:741-749 (2024). doi:10.1136/jmg-2024-109896

Rare-disease studies routinely begin with a number: the prevalence of the
disease. Duchenne affects roughly one in several thousand male births; other
disorders may affect one in tens or hundreds of thousands. After working with
these diseases, I became increasingly interested in something that seemed
almost too simple to ask.

What does prevalence mean in evolutionary terms? Why does one genetic disease
remain extremely rare while another complex disease affects millions of people?

## From disease genetics to population genetics

When I began working with Lin Chen and became more deeply involved in
statistical genetics and genomic epidemiology, I returned to this question.

It became apparent that the question could not be answered by experimental
manipulation of a model organism. What made it tractable was precisely the
feature of human genetics that is sometimes regarded as its weakness: we cannot
design the experiment.

Instead, millions of naturally occurring genetic perturbations have already
happened. Recombination has randomized them across generations. Natural
selection has acted on them. Physicians have recorded their phenotypic
consequences. Biobanks have measured disease and reproduction in hundreds of
thousands of individuals. GWAS and molecular QTL studies have connected those
phenotypes back to genetic variation.

The experiment, in a sense, has been running for generations. The challenge is
to extract its signal.

This led to a simple evolutionary argument. Natural selection does not act on
disease prevalence itself. Ultimately, it acts through fitness. A
disease-associated allele does not disappear merely because we call it a risk
allele. Its evolutionary trajectory depends on all of its phenotypic
consequences, including its effects on reproduction.

If genetic effects that increase disease liability are coupled to greater
reproductive success, selection on reproduction may oppose selection against
disease. This is antagonistic pleiotropy. If disease-increasing effects are
instead coupled to reduced reproductive success, the two selective consequences
act in the same direction.

That gives a testable prediction:

**Variation in disease prevalence may partly reflect how disease-associated
genetic variation is coupled to reproductive fitness.**

## Testing the idea

The first stage of this work characterized disease-fertility genetic coupling
across 97 complex diseases. The broader question became testable when
large-scale disease prevalence data could be integrated with genomic estimates
of reproductive effects.

Across thousands of diseases, a systematic pattern emerged. Disease prevalence
increased with more antagonistic disease-fertility coupling: genetic effects
associated with increased disease liability were also associated with greater
reproductive success. Conversely, prevalence decreased with synergistic
coupling, where disease-increasing effects were associated with reduced
reproductive success.

The same evolutionary logic also led us to examine differences between male and
female reproductive effects, revealing widespread sexually antagonistic effects
among disease-associated genes.

What began as a question raised by rare-disease genetics therefore became a
population-level evolutionary hypothesis: the reproductive consequences of
disease-associated genetic variation may help explain why some human diseases
are common while others remain rare.

This repository contains the code and analyses developed to test that
hypothesis.

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
