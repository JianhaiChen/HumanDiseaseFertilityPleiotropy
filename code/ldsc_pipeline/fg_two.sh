#!/bin/bash
cd /scratch/jianhaichen/selection
source /home/jianhaichen/.conda/etc/profile.d/conda.sh 2>/dev/null; conda activate ldsc 2>/dev/null||source activate ldsc
LD=eur_w_ld_chr/; HM=eur_w_ld_chr/w_hm3.snplist; M=ldsc/munge_sumstats.py; L=ldsc/ldsc.py
declare -A U=( [fgSepsis]='https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_AB1_OTHER_SEPSIS.gz' [fgVitD]='https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_E4_VIT_D_DEF.gz' )
declare -A NN=( [fgSepsis]=456181 [fgVitD]=467756 )
for lab in fgSepsis fgVitD; do
  wget -q -O $lab.gz "${U[$lab]}"
  zcat $lab.gz | awk -v N=${NN[$lab]} -f clean6_fix.awk > clean2/$lab.clean
  bash align_one.sh clean2/$lab.clean $HM al2/$lab.al.clean
  python $M --sumstats al2/$lab.al.clean --signed-sumstats BETA,0 --out munged2/$lab >/dev/null 2>&1
  echo "$lab aligned=$(wc -l<al2/$lab.al.clean) h2: $(python $L --h2 munged2/$lab.sumstats.gz --ref-ld-chr $LD --w-ld-chr $LD --out h2_$lab >/dev/null 2>&1; grep -i 'Total Observed' h2_$lab.log|head -1) | $(grep -i 'Mean Chi' h2_$lab.log|head -1)"
  for ref in father mother lifespan; do python $L --rg munged2/$lab.sumstats.gz,$ref.al.sumstats.gz --ref-ld-chr $LD --w-ld-chr $LD --out rglog3/rg_${lab}_$ref >/dev/null 2>&1; done
  echo "$lab father=$(grep -i 'Genetic Correlation:' rglog3/rg_${lab}_father.log|head -1|sed 's/.*: //') mother=$(grep -i 'Genetic Correlation:' rglog3/rg_${lab}_mother.log|head -1|sed 's/.*: //')"
done
echo FG_TWO_DONE
