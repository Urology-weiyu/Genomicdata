# ============================================================
# Title:    Venn Diagram of CUT&Tag Peaks and DEGs
# Purpose:  Identify overlap between CUT&Tag peaks and differential genes
# ============================================================

setwd("/path/to/workdir")

# ---------- [Read and filter DEGs] ----------
smchd1 <- read.csv("Limma_KOvsCon_SMCHD1_WithTPM.csv")
smchd1 <- subset(smchd1, adj.P.Val < 0.05 & logFC > 1)
write.csv(smchd1, "KOvsCon_SMCHD1_volcano_genes.csv", row.names = FALSE)

prmt5 <- read.csv("Limma_KOvsCon_PRMT5_WithTPM.csv")
prmt5 <- subset(prmt5, adj.P.Val < 0.05 & logFC > 1)
write.csv(prmt5, "KOvsCon_PRMT5_volcano_genes.csv", row.names = FALSE)

# ---------- [Read CUT&Tag peaks] ----------
bed <- read.csv("cuttag.csv")
bed <- subset(bed, gene_biotype == "protein_coding" &
                     annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)"))
bed <- unique(bed$SYMBOL)

# ---------- [Prepare gene lists] ----------
geneList <- list(
  "CUT&Tag"   = bed,
  "sgSMCHD1"  = smchd1$genesymbol,
  "sgPRMT5"   = prmt5$genesymbol
)

# ---------- [Intersection] ----------
intersectGenes <- Reduce(intersect, geneList)
print(intersectGenes)

# ---------- [Draw Venn diagram] ----------
library(VennDiagram)
venn_ploy <- venn.diagram(
  x = geneList,
  fill = c("#616fa9", "#da7733", "#3f7d7d"),
  print.mode = "raw",
  cat.cex = 1.8,
  cex = 1.6,
  cat.default.pos = "outer",
  scaled = FALSE,
  filename = NULL
)

grid.draw(venn_ploy)

# ---------- [Save to PDF] ----------
pdf("Venn.pdf", width = 4, height = 4)
grid.draw(venn_ploy)
dev.off()
