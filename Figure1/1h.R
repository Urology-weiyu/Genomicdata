# ============================================================
# Title:    GSEA Analysis and Visualization for KO vs Control
# System:   R 4.3.1
# Purpose:  Perform GSEA using logFC-ranked genes
# ============================================================


library(clusterProfiler)
library(pathview)
library(enrichplot)
library(msigdbr)
library(heatmap3)
library(GseaVis)
# Set clusterProfiler download option
R.utils::setOption("clusterProfiler.download.method", 'auto')

# ---------- [Prepare gene list for GSEA] ----------
markers_for_gsea <- read.csv('Limma_KOvsCon_PRMT5.csv')
gene_list <- markers_for_gsea$logFC
names(gene_list) <- markers_for_gsea$genesymbol
gene_list <- sort(gene_list, decreasing = TRUE)

# ---------- [Load custom gene sets for GSEA] ----------
mf_df <- read.csv('GSEAgenelist.csv')

# ---------- [Run GSEA] ----------
gsea_results <- GSEA(gene_list, TERM2GENE = mf_df)
gsea_results

# Save results
write_csv(as.data.frame(gsea_results), "gsea_results.csv")

# ---------- [Visualization] ----------
pdf('Cuproptosis_Promoting_GSEA.pdf', width = 3, height = 3)
gseaNb(object = gsea_results,
       geneSetID = 'Cuproptosis Promoting',
       addPval = TRUE,
       pvalX = 0.01,
       pvalY = 0.35,
       pCol = 'black',
       pvalSize = 4,
       pHjust = 0,
       htHeight = 0.5,
       ht.legend = FALSE,
       rmPrefix = FALSE)
dev.off()

# ---------- [Session Info Summary] ----------
# R version: 4.3.1 (2023-06-16)
# Platform:  x86_64-w64-mingw32/x64 (Windows 11)
# Key packages: clusterProfiler_4.10.0 pathview_1.42.0
# enrichplot_1.22.0 msigdbr_7.5.1 heatmap3_1.1.9
# GseaVis_0.1.0 
