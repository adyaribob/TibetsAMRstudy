library(cowplot)
library(dplyr)
library(stringr)
library(ggplot2)

### tables of heatmap and  barchart of random forest
dat.fig7b.heatmap.randomforest.contig.genus <- read.csv("dat.fig7b.heatmap.randomforest.contig.genus.csv")
dat.fig7b.barchart.randomforest.gene.contig.genus <- read.csv("dat.fig7b.barchart.randomforest.gene.contig.genus.csv")

(fig.7b.heatmap.randomforest.contig.genus <- ggplot(data = dat.fig7b.heatmap.randomforest.contig.genus[grepl(pattern = "IF", 
                                                                                               x = dat.fig7b.heatmap.randomforest.contig.genus$PredictedVariable),], 
                                            aes(x = factor(PredictedVariable), y = predictor, fill = IncMSE)) + geom_tile() + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
  theme_classic()+ 
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      label <- gsub("above", "<span style='color:#A52A2A;'>above</span>", label)
      label <- gsub("below", "<span style='color:#29A6A6;'>below</span>", label)
      label <- gsub("neutral", "<span style='color:grey;'>neutral</span>", label)
      label
    })
  }) +
  ggh4x::facet_grid2(cols = vars(media), scales = "free", space = "free") +
  geom_text(aes(x = PredictedVariable, y = predictor, label=sig), color="black", size=4) + 
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    #strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 13, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    #strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) )

(hm.clean <- fig.7b.heatmap.randomforest.contig.genus +
  theme(axis.title.y = element_blank(), axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none")
)
#### barchart #####
(dat.fig7b.barchart.randomforest.gene.contig.genus$media <- factor(dat.fig7b.barchart.randomforest.gene.contig.genus$media, 
                                                       levels = c("water", "sediment"))
)
(up <- ggplot(data = dat.fig7b.barchart.randomforest.gene.contig.genus, aes(x = PredictedVariable, y = r2)) + 
  geom_bar(stat = "identity", fill = "red") + ylab("") + 
  ggh4x::facet_grid2(cols = vars(media), scales = "free", space = "free") +
  scale_y_continuous(expand = c(0.01, 0),  limits = c(0,1))+theme_classic()+ 
  geom_text(aes(label = sig, y = r2 + 0.05), 
            angle = 360, hjust = 0.5, vjust = 0.5, size = 5, color = "black") +
  #geom_text(
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), #axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), axis.ticks.y= element_blank(),
        legend.position="none") 
)

(fig.7b <- cowplot::plot_grid( up,heatmap.randomforest.contig.genus, align = "v", nrow = 2, rel_heights = c(0.3, 1.5))
)
