library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)

my_colors <- c(
  "CAMRA"      = "red",
  "LDM-med"  = "blue",  
  "microHIMA"    = "grey60",  
  "MarZIC"     = "green",  
  "multimedia" = "yellow", 
  "CRAmed"     = "orange",  
  "CMM"        = "darkgreen",
  "PERMANOVA-med" = "purple",
  "MODIMA"     = "pink",
  "MedTest"    = "skyblue"
)

benchmark_methods <- c("HIMA", "LDM", "MarZIC", "multimedia", "CRAmed")

target_levels <- c(
  "Complete Null",
  "Exposure-only",
  "Outcome-only",
  "Disjoint (Balanced +/-)",
  "Disjoint (Dominant +)"
)

###### Fig 2 Empirical FDR of taxon-level mediation tests ######

plot_fig2_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>%       
  filter(num2 > 0) %>%            # mediation_signal
  filter(method %in% benchmark_methods) %>%
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

draw_fdr_barplot <- function(data, d_val) {
  sub_data <- data %>% filter(d == d_val)
  
  ggplot(sub_data, aes(x = num2, y = fdr_mean, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", size = 0.2, width = 0.7) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey", size = 0.8) +
    facet_grid(p_lab ~ n_lab) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 0.8)) +
    scale_fill_manual(values = my_colors) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),       
      legend.position = "bottom",         
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0, face = "bold", size = 16), 
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 14),     
      legend.title = element_text(size = 16)     
    ) +
    labs(
      title = paste0("Balanced (d = ", d_val, ")"),
      x = "Number of True Mediators",
      y = "Empirical FDR",
      fill = "Method"
    )
}

p1 <- draw_fdr_barplot(plot_fig2_data, 0.5) + 
  labs(title = "a. Balanced +/- ")

p2 <- draw_fdr_barplot(plot_fig2_data, 0.9) + 
  labs(title = "b. Dominant +")

final_fig2 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")
