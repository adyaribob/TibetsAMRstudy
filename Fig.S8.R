
dat.figS8.boxplot.abundance.occurrence.gene.contig <- read.csv("dat.figS8.boxplot.abundance.occurrence.gene.contig.csv")

ggplot(dat.figS8.boxplot.abundance.occurrence.gene.contig,  aes(x = x, y = Occurrences, fill = media)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  #scale_x_discrete(expand = c(0.15, 0.05)) +
  theme_classic()+
  ggh4x::facet_grid2(cols = vars(media))+
  labs(
    y = "Occurrences") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "Sewage" = "#c88b5cff", "None" = "grey80"))+
  xlab("")+ 
  
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = element_text(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))

