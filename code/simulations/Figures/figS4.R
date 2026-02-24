###### Fig S4 Empirical FDR of taxon-level mediation tests including CAMRA ######

plot_S4_data <- taxon_level_summary %>%
  filter(alpha == 0.05) %>%       
  filter(num2 > 0) %>%            # mediation_signal
  mutate(
    num2 = as.factor(num2),       
    n_lab = paste0("n = ", n),    
    p_lab = paste0("p = ", p)
  )

p1 <- draw_fdr_barplot(plot_S4_data, 0.5) + 
  labs(title = "a. Balanced +/-")

p2 <- draw_fdr_barplot(plot_S4_data, 0.9) + 
  labs(title = "b. Dominant +")

final_S4 <- (p1 /plot_spacer()/ p2) + 
  plot_layout(guides = "collect",
              heights = c(1,0.05,1)) & 
  theme(legend.position = "bottom")
