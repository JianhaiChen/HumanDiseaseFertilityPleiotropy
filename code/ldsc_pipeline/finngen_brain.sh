#!/bin/bash
cd /scratch/jianhaichen/selection
source /home/jianhaichen/.conda/etc/profile.d/conda.sh 2>/dev/null; conda activate ldsc 2>/dev/null||source activate ldsc
LD=eur_w_ld_chr/; HM=eur_w_ld_chr/w_hm3.snplist; M=ldsc/munge_sumstats.py; L=ldsc/ldsc.py
echo DOWNLOAD_START
wget -q -O fg_brain.gz 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_BRAIN_EXALLC.gz'
echo "downloaded $(ls -la fg_brain.gz|awk '{print $5}') bytes"
echo '=== FinnGen 列名 ==='; zcat fg_brain.gz | head -1
zcat fg_brain.gz | awk -v N=380643 -f clean6_fix.awk > clean2/fgBrain.clean
bash align_one.sh clean2/fgBrain.clean $HM al2/fgBrain.al.clean
echo "fgBrain rows=$(wc -l<al2/fgBrain.al.clean)"
python $M --sumstats al2/fgBrain.al.clean --signed-sumstats BETA,0 --out munged2/fgBrain >/dev/null 2>&1
python $L --h2 munged2/fgBrain.sumstats.gz --ref-ld-chr $LD --w-ld-chr $LD --out h2_fgBrain >/dev/null 2>&1
echo "fgBrain h2: $(grep -i 'Total Observed' h2_fgBrain.log|head -1) | $(grep -i 'Mean Chi' h2_fgBrain.log|head -1)"
if [ -s munged2/fgBrain.sumstats.gz ]; then for ref in father mother lifespan; do python $L --rg munged2/fgBrain.sumstats.gz,$ref.al.sumstats.gz --ref-ld-chr $LD --w-ld-chr $LD --out rglog3/rg_fgBrain_$ref >/dev/null 2>&1; done
  echo "fgBrain father=$(grep -i 'Genetic Correlation:' rglog3/rg_fgBrain_father.log|head -1|sed 's/.*: //') mother=$(grep -i 'Genetic Correlation:' rglog3/rg_fgBrain_mother.log|head -1|sed 's/.*: //')"; fi
echo FGBRAIN_DONE
