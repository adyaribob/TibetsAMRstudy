library(cowplot)
library(ggplot2)

### tables of heatmap and  barchart of random forest
dat.fig7a.heatmap.randomforest.gene.contig.modules <- read.csv("dat.fig7a.heatmap.randomforest.gene.contig.modules.csv")
dat.fig7a.barchart.randomforest.gene.contig.modules <- read.csv("dat.fig7a.barchart.randomforest.gene.contig.modules.csv")

dat.fig7a.heatmap.randomforest.gene.contig.modules$group <- 
  factor(dat.fig7a.heatmap.randomforest.gene.contig.modules$group, 
         levels = c("Sewage ARGs", "Sewage VFGs", "Sewage MGEs",
                    "Water ARGs", "Water VFGs", "Water MGEs",
                    "Sediment ARGs", "Sediment VFGs", "Sediment MGEs",
                    "ARG-MGE contigs", "ARG-VFG contigs", "ARG-MGE-VFG contigs"))

dat.fig7a.heatmap.randomforest.gene.contig.modules$module <- 
  factor(dat.fig7a.heatmap.randomforest.gene.contig.modules$module, 
         levels = c("Module 1", "Module 2", "Module 3", "Module 4", "Module 5", "Module 6",
                    "Module 7", "Module 8", "Module 9", "Module 10", "Module 11" ))
dat.fig7a.heatmap.randomforest.gene.contig.modules$media <- factor(dat.fig7a.heatmap.randomforest.gene.contig.modules$media, levels = c("Water", "Sediment"))


(fig.7a.random.forest.heatmap.gene.contig.module <- ggplot(data = dat.fig7a.heatmap.randomforest.gene.contig.modules, 
                                    aes(x = factor(group), y = module, fill = IncMSE)) + geom_tile() + 
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 15), legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
  geom_text(aes(x = group, y = module, label=sig), color="black", size=5) + 
  ggh4x::facet_grid2(cols = vars(media) ) + 
  theme_classic()+ 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.y.left  = element_text(size = 10, colour = "black"),
    #strip.text.x = element_text(size = 10.5, colour = "black"), 
    strip.text.y = element_text(size = 9, colour = "black"), 
    axis.title.y = element_text(size = 10),
    #strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(size = 10.5 ,face = "plain", hjust = 1, vjust = 1, angle = 45)))

##### random forest barchart #######
dat.fig7a.barchart.randomforest.gene.contig.modules$group <- 
  factor(dat.fig7a.barchart.randomforest.gene.contig.modules$group, 
         levels = c("Sewage ARGs", "Sewage VFGs", "Sewage MGEs",
                    "Water ARGs", "Water VFGs", "Water MGEs",
                    "Sediment ARGs", "Sediment VFGs", "Sediment MGEs",
                    "ARG-MGE contigs", "ARG-VFG contigs", "ARG-MGE-VFG contigs"))
dat.fig7a.barchart.randomforest.gene.contig.modules$media <- factor(dat.fig7a.barchart.randomforest.gene.contig.modules$media, levels = c("Water", "Sediment"))

(fig.7a.random.forest.barchart <- ggplot(data = dat.fig7a.barchart.randomforest.gene.contig.modules, aes(x = group, y = Rsquare)) + 
  geom_bar(stat = "identity", fill = "red") +
  ggh4x::facet_grid2(cols = vars(media) ) + 
  theme_classic()+ 
  geom_text(aes(x = group, y = Rsquare , label=sig), color="black", size=5) + 
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 7.5) +
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y.left  = element_text(size = 10, colour = "black"),
    strip.text.y = element_text(size = 9, colour = "black"), 
    axis.title.y = element_text(size = 10),
    legend.key = NULL,
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank() 
  )
)

cowplot::plot_grid(fig.7a.random.forest.barchart, fig.7a.random.forest.heatmap.gene.contig.module, 
          nrow = 2, cols = 1, align = "v", rel_widths = c(3,3), 
          rel_heights = c(1.2, 3))

