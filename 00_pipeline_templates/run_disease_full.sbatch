#!/bin/bash
# Stage 2, disease arm. One task per endpoint, all tissues in that task.
#SBATCH --job-name=dzfull
#SBATCH -p tier1q,tier2q,tier3q
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -o ${WORKDIR}/logs_dzfull/%A_%a.out
#SBATCH -e ${WORKDIR}/logs_dzfull/%A_%a.err
module load gcc 2>/dev/null; module load R/4.2.1 2>/dev/null
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BASEDIR}/prev_endpoints.tsv)
[ -z "$row" ] && exit 0
p=$(echo "$row" | cut -f1); f=$(echo "$row" | cut -f2)
cd ${WORKDIR}
export FUSIOS_PRECLUMP=preclump_full
export FUSIOS_OUTDIR=${WORKDIR}/tissue_sumtbl_full
Rscript code/fusios_iv5.R "$f" "$p" ""
