# ============================================================
# Title:    Volcano Plot for PRMT5 Knockdown RNA-seq
# System:   R 4.3.1
# Purpose:  Visualize differential expression results with
#           significance thresholds and highlight selected genes.
# ============================================================

# ---------- [Load Libraries & Data] ----------
library(ggplot2)
library(ggrepel)
df <- read.csv("Limma_KOvsCon_PRMT5.csv")

# ---------- [Define significance categories] ----------
df$change <- ifelse(df$adj.P.Val < 0.05 & df$logFC > 1, "sgPRMT5 enriched",
                    ifelse(df$adj.P.Val < 0.05 & df$logFC < -1, "sgPRMT5 depleted",
                           "Not significant"))

df$change <- factor(df$change, levels = c("sgPRMT5 enriched", "Not significant", "sgPRMT5 depleted"))

# ---------- [Select top genes for labeling] ----------
gene1 <- c("FDX1","PDHB","SLC25A3","SLC31A2","SLC31A1","DLAT","PDHA1","DLD","DLST",
           "LIAS","LIPT1","LIPT2","GCSH","MPC1")
gene2 <- c("ATP7A","ATP7B","MTF1","CDKN2A","GLS","MT2A")

top_up5_genes <- subset(df, genesymbol %in% c(gene1, gene2) & logFC > 1)

# ---------- [Volcano Plot] ----------
ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("sgPRMT5 enriched" = "#C6295C",
                                "sgPRMT5 depleted" = "#2C6DB2",
                                "Not significant" = "grey"),
                     guide = guide_legend(nrow = 3)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top_up5_genes, aes(label = genesymbol),
                  color = "black", size = 4, direction = "y", max.overlaps = Inf) +
  labs(
    x = expression(atop("sgPRMT5 versus sgNC", log[2](Fold~Change))),
    y = expression(-log[10]('adjusted P value')),
    color = ""
  ) +
  theme_test() +
  theme(
    axis.text = element_text(color = "black", size = 11),
    legend.position = c(0.25, 0.8),
    legend.text = element_text(size = 11)
  )

# ---------- [Save Plot] ----------
ggsave("sgPRMT5_volcano.pdf", width = 3, height = 3)

