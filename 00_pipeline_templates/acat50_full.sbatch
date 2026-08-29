#!/bin/bash
# Stage 3, disease arm.
#SBATCH -J mr_acat50
#SBATCH -p tier1q,tier2q,tier3q,express
#SBATCH --mem=24G
#SBATCH -t 05:00:00
#SBATCH --array=1-500%120
#SBATCH -o ${BASEDIR}/acat50_logs/%a.out
cd ${BASEDIR}
/apps/software/gcc-12.1.0/R/4.4.2/bin/Rscript acat50_full.R ${SLURM_ARRAY_TASK_ID} 500 2>&1 | grep -vE "built under|^Loading|package"
