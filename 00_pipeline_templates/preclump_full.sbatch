#!/bin/bash
# Stage 1. Array over the (tissue, chromosome) grid in preclump_jobs.tsv.
#SBATCH --job-name=pcfull
#SBATCH -p tier2q
#SBATCH --mem-per-cpu=40G
#SBATCH -t 16:00:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -o ${WORKDIR}/logs_pcfull/%A_%a.out
#SBATCH -e ${WORKDIR}/logs_pcfull/%A_%a.err
#SBATCH --open-mode=append
module load gcc 2>/dev/null; module load R/4.2.1 2>/dev/null
n=$((OFF + SLURM_ARRAY_TASK_ID))
row=$(sed -n "${n}p" ${WORKDIR}/preclump_jobs.tsv); [ -z "$row" ] && exit 0
t=$(echo "$row"|cut -f1); c=$(echo "$row"|cut -f2)
o=${WORKDIR}/preclump_full/$t/chr$c.txt
[ -s "$o" ] && exit 0
Rscript ${WORKDIR}/code/preclump_R_full.R "$t" "$c"
