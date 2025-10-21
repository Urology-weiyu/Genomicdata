# Drug combination viability heatmap generation
# Run on R 

setwd("/path/to/workdir")
rt <- read.csv("drug.csv")

# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(cowplot)

# Data reshape
rt_long <- rt %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Time",
    values_to = "Value"
  ) %>%
  mutate(
    Time = gsub("X", "", Time),
    Time = factor(Time, levels = rev(c("0","0.1","0.2","0.5","1","2"))),
    group = factor(group, levels = unique(group)),
    sample = factor(sample, levels = unique(sample))
  )

# Assign color to each inhibitor group
group_df <- rt_long %>%
  distinct(sample, group) %>%
  mutate(group_color = case_when(
    group == "DMSO"   ~ "#B0BEC5",
    group == "EZH2"   ~ "#f9766e",
    group == "DOT1L"  ~ "#fbbab6",
    group == "LSD1"   ~ "#e1c548",
    group == "GLP"    ~ "#f0e2a3",
    group == "KDM"    ~ "#5fa664",
    group == "HDAC"   ~ "#abd0a7",
    group == "P300"   ~ "#ca6a6b",
    group == "BRD"    ~ "#e5b5b5",
    group == "DNMT"   ~ "#4e79a6",
    group == "SMARCA" ~ "#bac4d0",
    group == "PRMT"   ~ "#45337f",
    TRUE ~ "grey"
  ))

rt_long <- rt_long %>%
  left_join(group_df, by = "sample")

# Define heatmap color palette
heat_colors <- colorRampPalette(c("#2166AC", "white", "#ffb30090", "#B2182B"))(100)

# Plot heatmap
p1 <- ggplot(rt_long, aes(sample, Time)) +
  annotate(
    "rect",
    xmin = as.numeric(factor(group_df$sample)) - 0.5,
    xmax = as.numeric(factor(group_df$sample)) + 0.5,
    ymin = -0.5,
    ymax = 0,
    fill = group_df$group_color
  ) +
  geom_tile(aes(fill = Value), color = "black") +
  scale_fill_gradientn(colors = heat_colors) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "Concentration (μM)", fill = "Cell Viability") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, size = 12, hjust = 1, vjust = .5, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 13, face = "bold")
  ) +
  coord_fixed() +
  guides(fill = guide_colorbar(
    position = "right",
    barwidth = .8,
    barheight = 8,
    direction = "vertical",
    ticks.colour = "black",
    ticks.outside = TRUE,
    frame.colour = "black"
  ))

ggsave("Fig1.png", plot = p1, width = 9, height = 4, dpi = 600)

# Plot group legend separately
group_color_df <- group_df %>%
  distinct(group, group_color)

legend_plot <- ggplot(group_color_df, aes(x = 1, y = group, fill = group)) +
  geom_point(shape = 21, size = 5) +
  scale_fill_manual(values = setNames(group_color_df$group_color, group_color_df$group)) +
  theme_void() +
  guides(fill = guide_legend(title = "Inhibitor", ncol = 2)) +
  theme(
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.8, "cm")
  )

legend <- cowplot::get_legend(legend_plot)
png("Fig1legend.png", width = 2, height = 3, units = "in", res = 300)
grid.newpage()
grid.draw(legend)
dev.off()


# ---------- [Session Info Summary] ----------
# R version: 4.3.1 (2023-06-16)
# Platform:  x86_64-w64-mingw32/x64 (Windows 11)
# Key packages: ggplot2_3.5.1, tidyr_1.3.0, dplyr_1.1.3,
#               viridis_0.6.4, RColorBrewer_1.1-3,
#               gridExtra_2.3, cowplot_1.1.1
