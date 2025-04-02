

dat.figS2.feast.barchart.water <- water

plot.feast.genes.water <- ggplot() + 
  geom_bar(data = dat.figS2.feast.barchart.water,
           mapping = aes(x = site, y = value, fill = variable), color = "black",
           stat = "identity",position = "stack")+
  theme_classic()+ ggh4x::facet_grid2(rows = vars(gene), cols = vars(river), 
                                      scales = "free", space = "free")+
  scale_y_continuous(
    limits = c(0, 1.05),          # Set y-axis range
    breaks = seq(0, 1, 0.2)    # Set labels at intervals of 0.2
  )+   
  labs(
    y = "FEAST source proportions") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "sewage" = "#c88b5cff", "unknown" = "grey80"), )+
  xlab("")+ 
  #stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
  #                   comparisons = pairs )+
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 0.5, angle = 90))


plot.feast.genes.sediment <- ggplot() + 
  geom_bar(data = dat.figS2.feast.barchart.sediment,
           mapping = aes(x = site, y = value, fill = variable), color = "black",
           stat = "identity",position = "stack")+
  theme_classic()+ ggh4x::facet_grid2(rows = vars(gene), cols = vars(river), 
                                      scales = "free", space = "free")+
  scale_y_continuous(
    limits = c(0, 1.05),          # Set y-axis range
    breaks = seq(0, 1, 0.2)    # Set labels at intervals of 0.2
  )+   
  labs(
    y = "FEAST source proportions") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "sewage" = "#c88b5cff", "unknown" = "grey80"), )+
  xlab("")+ 
  #stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
  #                   comparisons = pairs )+
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 0.5, angle = 90))



plot_grid(plot.feast.genes.water, plot.feast.genes.sediment,  align = "V", cols = 1,
          rel_widths = c(1, 1),
          rel_heights = c(1, 1)
