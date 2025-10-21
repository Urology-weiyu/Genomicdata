# ============================================================
# Title:    Volcano Plot for SMCHD1 Knockdown RNA-seq
# System:   R 4.3.1
# Purpose:  Visualize differential expression results with
#           significance thresholds and highlight selected genes.
# ============================================================

# ---------- [Load Libraries & Data] ----------
library(ggplot2)
library(ggrepel)
df <- read.csv("Limma_KOvsCon_SMCHD1.csv")

# ---------- [Define significance categories] ----------
df$change <- ifelse(df$adj.P.Val < 0.05 & df$logFC > 1, "sgSMCHD1 enriched",
                    ifelse(df$adj.P.Val < 0.05 & df$logFC < -1, "sgSMCHD1 depleted",
                           "Not significant"))

df$change <- factor(df$change,levels = c("sgSMCHD1 enriched","Not significant","sgSMCHD1 depleted"))

# ---------- [Select top genes for labeling] ----------
annogene <- subset(df,genesymbol %in% c('FDX1','SLC25A3','SLC31A1'))

# ---------- [Volcano Plot] ----------
ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),color = change)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("sgSMCHD1 enriched" = "#C6295C", "sgSMCHD1 depleted" = "#2C6DB2",
                                "Not significant" = "grey"),
                     guide = guide_legend(nrow = 3)) +
  scale_size_continuous(range = c(1,3)) +
  geom_vline(xintercept = c(-1,1),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",color = "grey50") +
  guides(size=guide_legend(order = 1)) +
  geom_text_repel(data = annogene, aes(label = genesymbol), 
                  color = 'black', size = 4, 
                  direction = "y", max.overlaps = Inf) +
  labs(x = expression(atop("sgSMCHD1 versus sgNC",log[2](Fold~Change)),
                      y = expression( -log[10]('adjusted P value'))),
       color = "")+
  theme_test() +
  theme(axis.text=element_text(color="black",size=11),
        legend.position = c(0.25,0.8),
        legend.text = element_text(size=11))

# ---------- [Save Plot] ----------
ggsave('sgSMCHD1_volcano.pdf',width = 3,height = 3) 
