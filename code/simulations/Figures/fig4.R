###### Fig4 FDR–power tradeoff for taxon-level mediator discovery ######

TARGET_FDR <- 0.05

fig4_data <- taxon_level_summary %>%
  filter(num2 > 0) %>% 
  group_by(template, n, p, d, num1A, num1B, num2, method) %>%
  summarise(
    valid_configs = sum(fdr_mean <= TARGET_FDR, na.rm = TRUE),
    max_power = if (valid_configs > 0) {
      max(power_mean[fdr_mean <= TARGET_FDR], na.rm = TRUE)
    } else {
      NA_real_ 
    },
    .groups = "drop"
  )

plot_fig4_data <- fig4_data %>%
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

draw_fig4_aligned <- function(data, d_val) {
  sub_data <- data %>% filter(d == d_val)
  ggplot(sub_data, aes(x = num2, y = max_power, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", size = 0.2, width = 0.7) +
    facet_grid(p_lab ~ n_lab) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1)) + 
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
      title = paste0("Max Power (FDR <= 0.05) | d = ", d_val),
      y = "Power@FDR<=0.05 ",
      x = "Number of True Mediators", 
      fill = "Method"
    )
}

p1 <- draw_fig4_aligned(plot_fig4_data, d_val = 0.5) + 
  labs(title = "a. Balanced +/-")

p2 <- draw_fig4_aligned(plot_fig4_data, d_val = 0.9) + 
  labs(title = "b. Dominant +")

final_fig4 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")
