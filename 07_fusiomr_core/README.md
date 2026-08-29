Model code sourced at run time by 01_transcriptome_wide_mr/run_fusiom_fast.R
(FusioMR_m) and 04_prevalence_finngen/fusios_iv5.R (FusioMR_s). Without these,
both drivers stop at their source()/sourceCpp() calls. Point ${FUSIOMR_CODE} here.

  FusioMR_s               FusioMR_m
    init_setup_seso.R       init_setup_semo.R
    gibbs_seso_uhp_only.cpp gibbs_semo_uhp_only_bound.cpp
    functions_perchr.R      set_variance_prior.R

The .cpp files are the Gibbs samplers, compiled at run time by Rcpp::sourceCpp().
sourceCpp() rewrites its cache index on every call, so a cacheDir shared across
concurrent array tasks corrupts it and the tasks die without a useful message.
Key the cache on $TMPDIR.

These routines implement the FusioMR method of Kang et al. (2026) and are included
so the analyses here run end to end. Cite the method's own release for the method.
