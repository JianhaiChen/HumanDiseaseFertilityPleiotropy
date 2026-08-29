#!/bin/bash
# Stage 2, chunked. Tissues split by chunks.txt; submit one array per CHUNK.
#SBATCH -p tier1q,tier2q,tier3q
#SBATCH --cpus-per-task=16
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -t 08:00:00
#SBATCH -o ${WORKDIR}/logs_dzchunk/%x_%A_%a.out
#SBATCH -e ${WORKDIR}/logs_dzchunk/%x_%A_%a.err
module load gcc 2>/dev/null; module load R/4.2.1 2>/dev/null
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BASEDIR}/prev_endpoints.tsv)
[ -z "$row" ] && exit 0
p=$(echo "$row" | cut -f1); f=$(echo "$row" | cut -f2)
TL=$(sed -n "${CHUNK}p" ${WORKDIR}/chunks.txt)
[ -z "$TL" ] && { echo "CHUNK=$CHUNK 取不到组织列表"; exit 1; }
cd ${WORKDIR}
export FUSIOS_PRECLUMP=preclump_full
export FUSIOS_OUTDIR=${WORKDIR}/tissue_sumtbl_full
Rscript code/fusios_iv5.R "$f" "$p" "$TL"
