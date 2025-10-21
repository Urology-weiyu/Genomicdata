# ============================================================
# Title:    Volcano Plot of Epigenetic Regulators in Adenoma/Polyp vs Mucosa
# Purpose:  Highlight Arginine/Histone methyltransferases, demethylases, and acetyltransferases
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ---------- [Read fold change data] ----------
tmp  <- fread('GSE207949_AdeCa-vs-Mucosa_foldchanges.csv')
tmp1 <- fread('GSE207949_Polyp-vs-Mucosa_foldchanges.csv')

# ---------- [Filter positive DEGs] ----------
tmp  <- subset(tmp, symbol != '' & padj != '' & log2FoldChange > 0)
tmp$Group <- 'tp1'
tmp1 <- subset(tmp1, symbol != '' & padj != '' & log2FoldChange > 0)
tmp1$Group <- 'tp2'

# ---------- [Combine datasets and calculate -log10(padj)] ----------
df2 <- rbind(tmp, tmp1)
df2$logp <- -log10(df2$padj)

# Flip log2FoldChange for tp2 to plot on left side
df2[grepl("tp2", Group), log2FoldChange := -log2FoldChange]

# ---------- [Load epigenetic regulator gene lists] ----------
gene <- fread('Epigenetic regulators.csv')
gene1 <- subset(gene, Function == 'Arginine Methyltransferase')$Gene
gene2 <- subset(gene, Function == 'Histone Lysine Methyltransferase')$Gene
gene3 <- subset(gene, Function == 'Lysine Demethylases')$Gene
gene4 <- subset(gene, Function == 'Histone Acetyltransferases')$Gene
gene5 <- subset(gene, Function == 'Histone Deacetylases')$Gene

# ---------- [Volcano plot] ----------
ggplot(df2, aes(x = log2FoldChange, y = logp)) +
  # Background points
  geom_point(data = subset(df2, !symbol %in% gene), shape = 21, color = "grey80", size = 2) +
  
  # Dashed lines for selected key genes
  geom_line(data = subset(df2, symbol %in% c('PRMT5','HDAC2','HDAC8','KDM5D','EHMT2')),
            color = "grey20", linetype = "dashed", size = 0.8) +
  
  # Highlight epigenetic regulator genes
  geom_point(data = subset(df2, symbol %in% gene1), color = "#c62828", size = 2.5) +
  geom_point(data = subset(df2, symbol %in% gene2), color = "#9c27b0", size = 2.5) +
  geom_point(data = subset(df2, symbol %in% gene3), color = "#283593", size = 2.5) +
  geom_point(data = subset(df2, symbol %in% gene4), color = "#00838f", size = 2.5) +
  geom_point(data = subset(df2, symbol %in% gene5), color = "#f9a825", size = 2.5) +
  
  # Label selected genes
  geom_text_repel(data = subset(df2, symbol %in% c('PRMT5')), 
                  aes(label = symbol), color = "#c62828", size = 5, fontface = "italic",
                  arrow = arrow(ends="first", length = unit(0.01, "npc")),
                  box.padding = 0.2, point.padding = 0.5, segment.color = 'black',
                  segment.size = 0.3, force = 1, max.iter = 3000, direction = "x") +
  geom_text_repel(data = subset(df2, symbol %in% c('HDAC2','HDAC8')), 
                  aes(label = symbol), color = "#f9a825", size = 5, fontface = "italic",
                  arrow = arrow(ends="first", length = unit(0.01, "npc")),
                  box.padding = 0.2, point.padding = 0.5, segment.color = 'black',
                  segment.size = 0.3, force = 1, max.iter = 3000, direction = "x") +
  geom_text_repel(data = subset(df2, symbol %in% c('EHMT2')), 
                  aes(label = symbol), color = "#9c27b0", size = 5, fontface = "italic",
                  arrow = arrow(ends="first", length = unit(0.01, "npc")),
                  box.padding = 0.2, point.padding = 0.5, segment.color = 'black',
                  segment.size = 0.3, force = 1, max.iter = 3000, direction = "x") +
  geom_text_repel(data = subset(df2, symbol %in% c('KDM5D')), 
                  aes(label = symbol), color = "#283593", size = 5, fontface = "italic",
                  arrow = arrow(ends="first", length = unit(0.01, "npc")),
                  box.padding = 0.2, point.padding = 0.5, segment.color = 'black',
                  segment.size = 0.3, force = 1, max.iter = 3000, direction = "x") +
  
  # Threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Labels and theme
  ylab(expression(-log[10](padj))) +
  xlab(expression(log[2]~FoldChange~"(versus Mucosa)")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 13, color = 'black'),
    axis.text.x = element_text(size = 13, color = 'black'),
    axis.title = element_text(size = 14),
    axis.ticks.length = unit(2, "mm")
  ) +
  scale_x_continuous(limits = c(-5.5, 5.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0, 0))
ggsave('FigS1b.pdf',width = 4,height = 4)
