Templates for the three compute stages, in run order. SLURM; the R scripts have
no scheduler dependency. Paths are variables: set ${BASEDIR} (working root) and
${WORKDIR} = ${BASEDIR}/fusios_finngen (instruments, per-tissue MR output).

Stage 1, instruments
  preclump_R_full.R      cis-eQTLs p<0.001, LD-clumped in R (100 kb, r2<0.1),
                         >=5 clumped instruments per gene
  preclump_full.sbatch   array over the (tissue, chromosome) grid

  No cap on eQTLs entering the clumping step. Capping first (top 30 by p) cuts
  surviving instruments several-fold and drops genes below the >=5 threshold.

Stage 2, transcriptome-wide MR
  run_disease_full.sbatch  one task per disease, all tissues
  run_dz_chunk.sbatch      tissues split into CHUNK groups, one array per chunk
                           (--export=ALL,CHUNK=k); shorter tasks, better packing
  run_fert_full.sbatch     fertility arm

  Both call 04_prevalence_finngen/fusios_iv5.R via FUSIOS_PRECLUMP (instruments)
  and FUSIOS_OUTDIR (output). Existing outputs are skipped, so a failed array can
  be resubmitted as is.

  Do not change --cpus-per-task mid-campaign: the RNG stream depends on mc.cores.
  Same for the reference panel.

  Sex-limited endpoints must not draw eQTLs from the other sex's tissue.
  fusios_iv5.R reads 05_sex_limited/endpoint_sex.tsv and skips ovary, uterus,
  vagina, fallopian tube and cervix for male-limited endpoints, testis and
  prostate for female-limited ones. Classify from the endpoint description in the
  FinnGen manifest, not the phenocode: codes are ambiguous for several anatomical
  terms (M13_CERVICALGIA is a neck; H8_EUSTSALP is an ear).

Stage 3, gene-level aggregation
  acat50_full.sbatch      disease arm
  acat_fert_full.sbatch   fertility arm

  Both call 04_prevalence_finngen/acat50.R. Signed gene-level effect = minimum-p
  tissue among those entering the ACAT. Cross-tissue IVW direction is reported as
  a sensitivity analysis.
