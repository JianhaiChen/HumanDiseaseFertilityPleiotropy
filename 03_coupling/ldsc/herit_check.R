suppressMessages(library(data.table))
w<-fread("why_ns_table.csv")
# 疾病遗传控制强度 (n_strong, sd_z1) 是否预测 |ρ_GE|?
cat("=== 疾病遗传控制 → |ρ_GE| 强度 (若疾病遗传度驱动,应正相关) ===\n")
cat(sprintf("  cor(n_strong, |rge|) = %.3f  P=%.3g\n",
  cor(w$n_strong_male,abs(w$rge_male)),cor.test(w$n_strong_male,abs(w$rge_male))$p.value))
cat(sprintf("  cor(sd_z1,    |rge|) = %.3f  P=%.3g\n",
  cor(w$sd_z1_male,abs(w$rge_male)),cor.test(w$sd_z1_male,abs(w$rge_male))$p.value))
# 分解非显著20个: 男女异号 vs 真近零
ns<-w[sig==F]
cat(sprintf("\n=== 非显著20个疾病的分解 ===\n"))
cat(sprintf("  男女异号(性别拮抗,均值抵消): %d\n",sum(ns$sex_discord)))
cat(sprintf("  男女同号但都弱(真·繁殖近正交): %d\n",sum(!ns$sex_discord)))
cat(sprintf("  这20个 |rge_avg| 中位 %.4f (确实平坦)\n",median(abs(ns$rge_avg))))
