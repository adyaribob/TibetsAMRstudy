library(ggplot2)


dat.figS3.barchart.abundance.occurrence.gene.contigs <- read.csv("dat.figS3.barchart.abundance.occurrence.gene.contigs.csv")

(barchart.percentage.neutral.contig.gene <- ggplot(data = dat.figS3.barchart.abundance.occurrence.gene.contigs,
                                                  aes(x = media, y = value, fill = neutral)) + 
  geom_bar(data = dat.figS3.barchart.abundance.occurrence.gene.contigs,
           mapping = aes(x = media, y = value, fill = neutral), 
           stat = "identity",position = "stack")+
  theme_classic()+ ggh4x::facet_grid2(~level , 
                                      scales = "free", space = "free",) + 
  scale_y_continuous(
    limits = c(0, 101),          # Set y-axis range
    breaks = seq(0, 101, 20),
    expand = c(0.01, 0)# Set labels at intervals of 0.2
  )+   
  labs(
    y = "Percentage (%)") +
  geom_text(data = dat.figS3.barchart.abundance.occurrence.gene.contigs, 
            aes(label = value), position = position_stack(vjust = 0.5)) +  # Add text labels
  #  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
  #                               "sewage" = "#c88b5cff", "none" = "grey80"), )+
  xlab("")+ 
  scale_fill_manual(name = "lefse", values = c("above" = "#A52A2A",  "neutral" = "grey",
                                               "below" = "#29A6A6" ))+
  #stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
  #                   comparisons = pairs )+
  theme(ggh4x.facet.nestline = element_line(linetype = 3),
        axis.text.y.left  = element_text(size = 12, colour = "black"),
        strip.text.x = element_text(size = 12, colour = "black"), 
        strip.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y = element_text(size = 12),
        strip.background = element_rect(size = 0.5),
        legend.key = NULL,
        legend.position = "bottom",
        axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))
)
