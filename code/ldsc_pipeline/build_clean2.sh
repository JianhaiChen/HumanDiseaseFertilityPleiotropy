#!/bin/bash
F="$1"; N="$2"; OUT="$3"
zcat -f "$F" | awk -v N=$N 'BEGIN{OFS="\t"}
  function pick(c,  k,a,n){n=split(c,a,","); for(k=1;k<=n;k++) if(a[k] in h) return h[a[k]]; return 0}
  {gsub(/\r/,"")}
  NR==1{ for(i=1;i<=NF;i++) h[$i]=i
    s=pick("rsid,hm_rsid,rs_id,SNP,snp,MarkerName,markername,variant_id")
    e=pick("effect_allele,hm_effect_allele,A1,a1,EA,Allele1,ALT")
    o=pick("other_allele,non_effect_allele,hm_other_allele,A2,a2,a0,NEA,Allele2,REF")
    b=pick("beta,hm_beta,BETA,Beta,effect_size,Effect,EFFECT")
    orr=pick("odds_ratio,OR,hm_odds_ratio,OddsRatio")
    p=pick("p_value,hm_p_value,P,pvalue,p,pval,Pval,P_value,P_BOLT_LMM")
    print "SNP","A1","A2","BETA","P","N"; next }
  { rs=$s; if(rs !~ /^rs/) next
    if(b!=0){ bv=$b } else if(orr!=0 && $orr+0>0){ bv=log($orr) } else next
    if(bv=="" || bv=="NA") next
    pv=$p; if(pv=="" || pv=="NA" || pv+0<=0) next
    print rs, toupper($e), toupper($o), bv, pv, N }' > "$OUT"
