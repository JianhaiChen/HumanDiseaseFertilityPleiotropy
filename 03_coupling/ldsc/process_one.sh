#!/bin/bash
dz="$1"; file="$2"; N="$3"
WD=${BASEDIR}/selection
GD=${LABDIR}/${USER}/gwasfullsummary_syn
HM=$WD/eur_w_ld_chr/w_hm3.snplist; M=$WD/ldsc/munge_sumstats.py
lab=$(echo "$dz" | tr -c 'A-Za-z0-9' '_')
[ -e "$GD/$file" ] || { echo "$lab NOFILE"; exit; }
bash $WD/build_clean2.sh "$GD/$file" "$N" "$WD/clean2/$lab.clean"
bash $WD/align_one.sh "$WD/clean2/$lab.clean" "$HM" "$WD/al2/$lab.al.clean"
nr=$(wc -l < "$WD/al2/$lab.al.clean" 2>/dev/null)
[ "${nr:-0}" -gt 100000 ] && python "$M" --sumstats "$WD/al2/$lab.al.clean" --signed-sumstats BETA,0 --out "$WD/munged2/$lab" >/dev/null 2>&1
[ -s "$WD/munged2/$lab.sumstats.gz" ] && echo "$lab OK rows=$nr" || echo "$lab FAIL rows=$nr"
