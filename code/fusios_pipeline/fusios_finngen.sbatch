#!/bin/bash
#SBATCH --job-name=fsFG
#SBATCH --partition=tier2q
#SBATCH --mem-per-cpu=25G
#SBATCH --time=6:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --error=logs/%A_%a.err
#SBATCH --output=logs/%A_%a.out
set -uo pipefail
cd /scratch/jianhaichen/fusios_finngen
module load gcc 2>/dev/null; module load R/4.2.1 2>/dev/null
LIST=$1
ph=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -z "$ph" ] && { echo "no phenocode at row $SLURM_ARRAY_TASK_ID"; exit 0; }
f=/scratch/jianhaichen/finngen_r13/finngen_R13_${ph}.gz
[ -s "$f" ] || { echo "missing $f"; exit 0; }
echo "TASK $SLURM_ARRAY_TASK_ID $ph"
Rscript code/fusios_finngen.R "$f" "$ph"
