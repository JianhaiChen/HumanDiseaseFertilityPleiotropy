suppressMessages({library(openxlsx); library(data.table)})
setwd("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait")
wb <- loadWorkbook("Supplementary Table.xlsx")
d  <- setDT(read.xlsx(wb, sheet="Supplementary Table 10", startRow=2))
nw <- fread("tableS10_fdr05_block.tsv")
drop <- intersect(setdiff(names(nw), "phenocode"), names(d))
if (length(drop)) d[, (drop) := NULL]
d <- merge(d, nw, by="phenocode", all.x=TRUE, sort=FALSE)
hdr <- read.xlsx(wb, sheet="Supplementary Table 10", rows=1, colNames=FALSE)[1,1]
desc <- paste0(trimws(hdr),
  " Columns ending _fdr05 recompute the directional-class counts and antagonistic",
  " shares at Benjamini-Hochberg FDR < 0.05 applied within trait to the 50-tissue",
  " gene-level ACAT P values; these are the values used in the main text and in",
  " Fig. 4a,b and Supplementary Figs. 10 and 13. Columns without the suffix retain",
  " the nominal P < 0.05 definition used in earlier versions.")
removeWorksheet(wb, "Supplementary Table 10")
addWorksheet(wb, "Supplementary Table 10")
writeData(wb, "Supplementary Table 10", desc, startRow = 1, colNames = FALSE)
writeData(wb, "Supplementary Table 10", d, startRow = 2)
ord <- suppressWarnings(as.numeric(gsub("[^0-9]", "", names(wb))))
worksheetOrder(wb) <- order(ord, na.last = FALSE)
saveWorkbook(wb, "Supplementary Table.xlsx", overwrite = TRUE)
cat("表10 已写回:", nrow(d), "行 x", ncol(d), "列\n")
cat("工作表:", paste(names(wb), collapse = " | "), "\n")
