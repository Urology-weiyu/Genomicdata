# ============================================================
# Title:    Correlation Plots
# ============================================================
setwd('path/to/workdir')
library(ggplot2)
library(ggpointdensity)
library(viridis)


format_p_value <- function(p){
  if(p < 0.001){
    return("P < 0.001")
  } else if(p < 0.01){
    return("P < 0.01")
  } else if(p < 0.05){
    return("P < 0.05")
  } else{
    return(paste0("P = ", signif(p, 3)))
  }
}
# ---------- [Spearman correlation] ----------
plot_correlation <- function(file, xlab, ylab, outpdf){
  tmp <- read.csv(file)
  cor_result <- cor.test(tmp[[1]], tmp[[2]], method = "spearman")
  r_value <- round(cor_result$estimate, 3)
  p_value <- cor_result$p.value
  label_text <- paste0("R = ", r_value, "\n", format_p_value(p_value))
  
  ggplot(tmp, aes(x = tmp[[1]], y = tmp[[2]])) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    labs(x = xlab, y = ylab) +
    annotate("text", x = 150, y = max(tmp[[2]]), 
             label = label_text, hjust = 0, vjust = 1, size = 6) +
    theme_classic(base_size = 15) +
    theme(axis.text = element_text(color = "black", size = 14),
          legend.position = "none",
          panel.grid = element_blank())
  
  ggsave(outpdf, width = 3.5, height = 3.5)
}

plot_correlation('Figure3g.csv', 'H score of SMCHD1', 'H score of FDX1', 'Fig3g.pdf')
plot_correlation('Figure3j.csv', 'H score of SMCHD1', 'H score of SLC31A1', 'Fig3j.pdf')
plot_correlation('Figure3m.csv', 'H score of SMCHD1', 'H score of SLC25A3', 'Fig3m.pdf')
