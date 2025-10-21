# ============================================================
# Title:    ssGSEA Analysis of Cuproptosis in FAP vs CRC
# Purpose:  Calculate Cuproptosis ssGSEA scores and compare between FAP and CRC
# ============================================================

setwd('/path/to/workdir')

# ---------- [Load libraries] ----------
library(data.table)
library(dplyr)
library(GSEABase)
library(clusterProfiler)
library(GSVA)
library(ggpubr)
library(ggplot2)

# ---------- [Read expression data] ----------
exp <- fread('GSE207949_norm_counts_TPM_GRCh38.p13_NCBI.tsv')
exp <- aggregate(. ~ Symbol, max, data = exp)
rownames(exp) <- exp$Symbol
exp <- exp[, -1]

# ---------- [ssGSEA calculation] ----------
exp <- log2(exp + 0.1)
geneset <- read.gmt('cuproptosis_extend_GSEA.gmt')

es.max <- gsva(as.matrix(exp), geneset, method = 'ssgsea',
               kcdf = "Gaussian",
               mx.diff = FALSE, verbose = FALSE, 
               parallel.sz = 1)

group <- es.max %>% t() %>% as.data.frame()
group$id <- rownames(group)
colnames(group) <- c('Cupro', 'x')



# ---------- [Subset samples by Disease_State] ----------
cli=fread('cli.csv')
cli=cli[,c(1,2,12)]
colnames(cli)=c('x','Patient','Disease_State')
tmp=merge(group,cli,by='x')
tmp$Disease_State=gsub('Stage: ','',tmp$Disease_State)

tmp <- group
crc <- subset(tmp, Disease_State == 'CRC')
fap <- subset(tmp, Disease_State == 'FAP')

# ---------- [Wilcoxon test] ----------
pvalue <- wilcox.test(crc$Cupro, fap$Cupro)$p.value

format_p_value <- function(p) {
  if (p < 0.001) {
    return("P < 0.001")
  } else if (p < 0.01) {
    return("P < 0.01")
  } else if (p < 0.05) {
    return("P < 0.05")
  } else {
    return(paste("P =", format(p, digits = 2)))
  }
}
formatted_p_values <- sapply(pvalue, format_p_value)

# ---------- [Prepare plotting data] ----------
tmp <- subset(tmp, Disease_State %in% c('FAP','CRC'))
tmp$Disease_State <- factor(tmp$Disease_State, levels = c('FAP','CRC'))

# ---------- [Plot boxplot with significance] ----------
ggplot(tmp, aes(x = Disease_State, y = Cupro, fill = Disease_State)) +
  stat_boxplot(geom = 'errorbar', size = 0.5, width = 0.4) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  theme_classic() +
  labs(y = 'Cuproptosis score', x = '') +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 13, color = 'black'),
    axis.text.x = element_text(size = 13, color = 'black', angle = 45, vjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.ticks.length = unit(2, "mm"),
    legend.position = 'none'
  ) +
  scale_y_continuous(limits = c(8.5, 10), expand = c(0, 0)) +
  scale_fill_manual(values = c('FAP' = '#ffcc80', 'CRC' = '#e57373')) +
  geom_signif(
    comparisons = list(c('CRC', 'FAP')),
    annotations = formatted_p_values,
    tip_length = 0,
    vjust = -0.5,
    step_increase = 0.15
  )
ggsave('FigS1a.pdf',width = 2,height = 3.4)

# ---------- [Session Info Summary] ----------
# R version: 4.3.1 (2023-06-16)
# Platform:  x86_64-w64-mingw32/x64 (Windows 11)
# Key packages: clusterProfiler_4.10.0 GseaVis_0.1.0
# clusterProfiler_4.10.0 GSVA_1.50.0
