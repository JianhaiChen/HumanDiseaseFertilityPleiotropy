#!/bin/bash
awk 'BEGIN{OFS="\t"}
  NR==FNR{a1[$1]=$2;a2[$1]=$3;next}
  FNR==1{print "SNP","A1","A2","BETA","P","N";next}
  ($1 in a1){
    if($2==a1[$1]&&$3==a2[$1]) print $1,$2,$3,$4,$5,$6
    else if($2==a2[$1]&&$3==a1[$1]) print $1,a1[$1],a2[$1],(0-$4),$5,$6
  }' "$2" "$1" > "$3"
