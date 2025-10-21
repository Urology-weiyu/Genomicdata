# ============================================================
# Title:    RNA-seq Differential Expression Analysis
# System:   R 4.3.1
# Purpose:  Perform differential expression analysis using
#           easyTCGA/limma on raw counts.
# ============================================================

# ---------- [Load Libraries & Data] ----------
library(easyTCGA)
count <- read.csv("count_file.csv")
rownames(count) <- count$Gene_id
count <- count[ , -1]   # remove Gene_id column
condition <- factor(c(rep("Con",3), rep("KO",3)), levels = c("Con","KO"))

# ---------- [Differential Expression Analysis] ----------
diff_res <- diff_analysis(
  exprset = count,
  group = condition,
  is_count = TRUE
)

res <- diff_res$deg_limma
# ---------- [Export Results] ----------
write.csv(res, "Limma_KOvsCon.csv", row.names = FALSE)
