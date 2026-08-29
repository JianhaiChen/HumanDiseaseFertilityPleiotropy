#!/bin/bash
# Stage 3, fertility arm.
#SBATCH -J acatfertfull
#SBATCH -p tier2q
#SBATCH --mem=24G
#SBATCH -t 06:00:00
#SBATCH --array=1-2
#SBATCH -o ${BASEDIR}/acat_fert_full_logs/%a.out
cd ${BASEDIR}
/apps/software/gcc-12.1.0/R/4.4.2/bin/Rscript acat50_fert_full.R ${SLURM_ARRAY_TASK_ID} 2
