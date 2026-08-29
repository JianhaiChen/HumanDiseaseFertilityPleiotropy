#!/bin/bash
cd ${BASEDIR}/selection
source /home/${USER}/.conda/etc/profile.d/conda.sh 2>/dev/null; conda activate ldsc 2>/dev/null || source activate ldsc
GD=${LABDIR}/${USER}/gwasfullsummary_syn
HM=eur_w_ld_chr/w_hm3.snplist; LD=eur_w_ld_chr/; M=ldsc/munge_sumstats.py; L=ldsc/ldsc.py
labs=(Cataract Melanoma Pancreatic_Cancer Panic non_Hodgkins_lymphoma thyroidcancer)
files=(Cataract.gwas_summary Melanoma.gwas_summary Pancreatic_Cancer.gwas_summary Panic.gwas_summary non-Hodgkins_lymphoma.gwas_summary Thyroid_Cancer.gwas_summary)
Ns=(586243 417127 411013 10240 412750 1620354)
for i in ${!labs[@]}; do
  lab=${labs[$i]}
  zcat -f "$GD/${files[$i]}" | awk -v N=${Ns[$i]} -f clean6_fix.awk > clean2/$lab.clean
  bash align_one.sh clean2/$lab.clean $HM al2/$lab.al.clean
  nr=$(wc -l < al2/$lab.al.clean)
  python $M --sumstats al2/$lab.al.clean --signed-sumstats BETA,0 --out munged2/$lab >/dev/null 2>&1
  ok=$([ -s munged2/$lab.sumstats.gz ] && echo MUNGED || echo FAIL)
  echo "$lab rows=$nr $ok"
  [ "$ok" = MUNGED ] || continue
  for ref in father mother lifespan; do python $L --rg munged2/$lab.sumstats.gz,$ref.al.sumstats.gz --ref-ld-chr $LD --w-ld-chr $LD --out rglog3/rg_${lab}_${ref} >/dev/null 2>&1; done
  echo "$lab father=$(grep -i 'Genetic Correlation:' rglog3/rg_${lab}_father.log|head -1|sed 's/.*: //') mother=$(grep -i 'Genetic Correlation:' rglog3/rg_${lab}_mother.log|head -1|sed 's/.*: //')"
done
echo ALL6C_DONE
