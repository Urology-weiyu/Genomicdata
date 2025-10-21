# ============================================================
# Title:    GO Enrichment Analysis of Upregulated Genes in KO vs Control
# Purpose:  Perform GO enrichment for significantly upregulated DEGs
# ============================================================

setwd('path/to/workdir')

# ---------- [Read DEG table and filter significant genes] ----------
tmp <- read.csv('Limma_KOvsCon_SMCHD1.csv')
tmp <- subset(tmp, P.Value < 0.05 & logFC > 0.1)

# ---------- [Load libraries] ----------
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tibble)
library(ggplot2)

# ---------- [Convert gene symbols to Entrez IDs] ----------
genelist <- bitr(
  tmp$genesymbol,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
tmp <- inner_join(tmp, genelist, by = c("genesymbol" = "SYMBOL"))

# ---------- [GO enrichment analysis] ----------
go <- enrichGO(
  gene          = tmp$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",          
  pAdjustMethod = "BH",          
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE          
)

# ---------- [Save results] ----------
go_res <- go@result
write.csv(go_res, file = "GO_SMCHD1.csv")


# ============================================================
# Title:    GO Enrichment Visualization (BP and CC)
# Purpose:  Barplot visualization of enriched GO terms
# ============================================================

# ---------- [Read GO result] ----------
df <- read.csv("pathways.csv", header = TRUE)
df$logP <- -log10(df$qvalue)
df <- df %>% arrange(desc(logP))
df$Description <- factor(df$Description, levels = rev(df$Description))

# ---------- [Filter BP terms and get top 5 genes] ----------
df_bp <- subset(df, ONTOLOGY == 'BP')
df_bp$top5_genes <- sapply(strsplit(df_bp$geneID, "/"),
                            function(x) paste(head(x, 5), collapse = "/"))

x_max <- max(df_bp$logP)

# ---------- [Barplot with circles indicating gene counts] ----------
p <- ggplot() +
  geom_bar(data = df_bp,
           aes(x = logP, y = Description, fill = ONTOLOGY),
           width = 0.75,
           stat = 'identity',
           alpha = 0.6) +
  geom_point(data = df_bp,
             aes(x = x_max + 1, y = Description),
             shape = 21, size = 7.5, fill = "#e57373", alpha = 0.5, color = 'white') +
  geom_text(data = df_bp,
            aes(x = x_max + 1, y = Description, label = Count),
            size = 3.5, color = "black") +
  scale_fill_manual(values = c('#e57373')) +
  scale_x_continuous(expand = c(0, 1)) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        legend.position = 'none') +
  geom_text(data = df_bp,
            aes(x = 0.3, y = Description, label = Description),
            size = 4.5, hjust = 0) +
  labs(x = expression(-log[10](qvalue)), y = "GO:BP Terms")
p
ggsave('GO_BP.pdf',height = 3.5,width = 5.5)


# ---------- [Filter CC terms and get top 5 genes] ----------
df_cc <- subset(df, ONTOLOGY == 'CC')
df_cc$top5_genes <- sapply(strsplit(df_cc$geneID, "/"),
                           function(x) paste(head(x, 5), collapse = "/"))

x_max <- max(df_cc$logP)

# ---------- [Barplot with circles indicating gene counts] ----------
p <- ggplot() +
  geom_bar(data = df_cc,
           aes(x = logP, y = Description, fill = ONTOLOGY),
           width = 0.75,
           stat = 'identity',
           alpha = 0.6) +
  geom_point(data = df_cc,
             aes(x = x_max + 2, y = Description),
             shape = 21, size = 7.5, fill = "#81c784", alpha = 0.5, color = 'white') +
  geom_text(data = df_cc,
            aes(x = x_max + 2, y = Description, label = Count),
            size = 3.5, color = "black") +
  scale_fill_manual(values = c('#81c784')) +
  scale_x_continuous(expand = c(0, 2)) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        legend.position = 'none') +
  geom_text(data = df_cc,
            aes(x = 0.4, y = Description, label = Description),
            size = 4.5, hjust = 0) +
  labs(x = expression(-log[10](qvalue)), y = "GO:CC Terms")
p
ggsave('GO_CC.pdf', height = 3.5, width = 5.5)
