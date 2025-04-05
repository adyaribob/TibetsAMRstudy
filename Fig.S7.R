
library(ggplot2)
### a ####
dat.figS7.barchart.contigs.number.based.on.neutral.lefse.pathogenicity.water <- read.csv("dat.figS7.barchart.contigs.number.based.on.neutral.lefse.pathogenicity.water.csv")
dat.figS7.barchart.euk.per.module.water <- read.csv("dat.figS7.barchart.euk.per.module.water.csv")
dat.figS7.heatmap.cor.subdivision.contigs.water <- read.csv("dat.figS7.heatmap.cor.subdivision.contigs.water.csv")

barchart.up <- dat.figS7.barchart.contigs.number.based.on.neutral.lefse.pathogenicity.water %>% 
  group_by(split) %>% summarise(n = n())
barchart.up$split <- gsub("NA_", replacement = "",x = barchart.up$split)
barchart.up$split <- gsub("non_pathogen", replacement = "non pathogen",x = barchart.up$split)
barchart.up$split <- gsub("_", replacement = " ",x = barchart.up$split)
barchart.up$split<-factor(barchart.up$split,levels=c("IF above non pathogen","IF above pathogen",
                                                     "IF below non pathogen","IF below pathogen",
                                                     "IF neutral non pathogen","IF neutral pathogen",
                                                     "WA above non pathogen","WA above pathogen",
                                                     "WA below non pathogen","WA below pathogen",
                                                     "WA neutral non pathogen","WA neutral pathogen",
                                                     "SD above non pathogen","SD above pathogen",
                                                     "SD below non pathogen","SD below pathogen",
                                                     "SD neutral non pathogen","SD neutral pathogen",
                                                     "above non pathogen","above pathogen",
                                                     "below non pathogen","below pathogen",
                                                     "neutral non pathogen","neutral pathogen"
))


barchart.right <- dat.figS7.barchart.euk.per.module.water %>% group_by(module.from) %>% summarise(n = n())
barchart.right$module.from <- factor(barchart.right$module.from, 
                                     levels = c("M_1", "M_2", "M_3", "M_4",
                                                "M_5", "M_6", "M_7", "M_8",
                                                "M_9", "M_10", "M_11"))


plot.heatmap.module.cor.contigs.lefse.neutral.water <- ggplot(dat.figS7.heatmap.cor.subdivision.contigs.water[dat.figS7.heatmap.cor.subdivision.contigs.water$group %in% barchart.up$split ,], 
                                                              aes(x = group, y = split)) + 
  geom_tile(aes(x = group, y = split, alpha = (mean),
                fill = (mean)))+ 
  theme_classic()+ 
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      label <- gsub("above", "<span style='color:#A52A2A;'>above</span>", label)
      label <- gsub("below", "<span style='color:#29A6A6;'>below</span>", label)
      label <- gsub("neutral", "<span style='color:grey;'>neutral</span>", label)
      label
    })
  }) +
  #ggh4x::facet_grid2(cols = vars(module.to))+
  # scale_alpha_manual(values = c(0.5, 0.75, 1)) +
  scale_fill_gradient2(low = "grey100",  mid = "grey50", high = "red", midpoint = 0.2) +
  theme_classic()+
  geom_text(
    aes(x = group, y = split, label=letters), 
    color="black", size = 4.5) + 
  #facet_grid(~facet, scales = "free", space = "free")+
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 10),
        axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 11, angle = 45),
        axis.text.y = ggtext::element_markdown(hjust = 1, vjust = 0.5, size = 11),
        strip.text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "bottom",
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("")+ylab("")


hm.clean <- plot.heatmap.module.cor.contigs.lefse.neutral.water +
  theme(axis.title.y = element_blank(), axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none")

right <- ggplot(data = barchart.right, aes(x = module.from, y = n)) + 
  geom_bar(stat = "identity") + ylab("") + 
  scale_y_continuous(expand = c(0.01, 0), limits = c(0,160))+theme_classic()+
  geom_text(aes(label = n, y = n + 5), 
            angle = 0, hjust = 0.5, vjust = 0.5, size = 5, color = "red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y= element_blank(),
        legend.position="none")  +coord_flip()
right

up <- ggplot(data = barchart.up, aes(x = split, y = n)) + 
  geom_bar(stat = "identity") + ylab("") + 
  scale_y_continuous(expand = c(0.01, 0),  limits = c(0,68))+theme_classic()+ 
  geom_text(aes(label = n, y = n + 2), 
            angle = 360, hjust = 0.5, vjust = 0.5, size = 5, color = "red") +
  #geom_text(
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y= element_blank(),
        legend.position="none") 
up

plot_grid(
  up, ggplot(), 
  hm.clean, right, 
  byrow = TRUE,
  nrow = 2, ncol = 2, rel_widths = c(1.5 ,0.3), rel_heights = c(0.3, 1.5)
)

### b ####
dat.figS7.barchart.contigs.number.based.on.neutral.lefse.pathogenicity.sediment <- read.csv("dat.figS7.barchart.contigs.number.based.on.neutral.lefse.pathogenicity.sediment.csv")
dat.figS7.barchart.euk.per.module.sediment <- read.csv("dat.figS7.barchart.euk.per.module.sediment.csv")
dat.figS7.heatmap.cor.subdivision.contigs.sediment <- read.csv("dat.figS7.heatmap.cor.subdivision.contigs.sediment.csv")


barchart.up <- dat.figS7.barchart.contigs.number.based.on.neutral.lefse.pathogenicity %>% 
  group_by(split) %>% summarise(n = n())
barchart.up$split <- gsub("NA_", replacement = "",x = barchart.up$split)
barchart.up$split <- gsub("non_pathogen", replacement = "non pathogen",x = barchart.up$split)
barchart.up$split <- gsub("_", replacement = " ",x = barchart.up$split)
barchart.up$split<-factor(barchart.up$split,levels=c("IF above non pathogen","IF above pathogen",
                                                     "IF below non pathogen","IF below pathogen",
                                                     "IF neutral non pathogen","IF neutral pathogen",
                                                     "WA above non pathogen","WA above pathogen",
                                                     "WA below non pathogen","WA below pathogen",
                                                     "WA neutral non pathogen","WA neutral pathogen",
                                                     "SD above non pathogen","SD above pathogen",
                                                     "SD below non pathogen","SD below pathogen",
                                                     "SD neutral non pathogen","SD neutral pathogen",
                                                     "above non pathogen","above pathogen",
                                                     "below non pathogen","below pathogen",
                                                     "neutral non pathogen","neutral pathogen"
))


barchart.right <- dat.figS7.barchart.euk.per.module.sediment %>% group_by(module.from) %>% summarise(n = n())
barchart.right$module.from <- factor(barchart.right$module.from, 
                                     levels = c("M_1", "M_2", "M_3", "M_4",
                                                "M_5", "M_6", "M_7", "M_8",
                                                "M_9", "M_10", "M_11"))


plot.heatmap.module.cor.contigs.lefse.neutral.sediment <- ggplot(dat.figS7.heatmap.cor.subdivision.contigs.sediment[dat.figS7.heatmap.cor.subdivision.contigs.sediment$group %in% barchart.up$split ,], 
                                                                 aes(x = group, y = split)) + 
  geom_tile(aes(x = group, y = split, alpha = (mean),
                fill = (mean)))+ 
  theme_classic()+ 
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      label <- gsub("above", "<span style='color:#A52A2A;'>above</span>", label)
      label <- gsub("below", "<span style='color:#29A6A6;'>below</span>", label)
      label <- gsub("neutral", "<span style='color:grey;'>neutral</span>", label)
      label
    })
  }) +
  #ggh4x::facet_grid2(cols = vars(module.to))+
  # scale_alpha_manual(values = c(0.5, 0.75, 1)) +
  scale_fill_gradient2(low = "grey100",  mid = "grey50", high = "red", midpoint = 0.2) +
  theme_classic()+
  geom_text(
    aes(x = group, y = split, label=letters), 
    color="black", size = 4.5) + 
  #facet_grid(~facet, scales = "free", space = "free")+
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 10),
        axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 11, angle = 45),
        axis.text.y = ggtext::element_markdown(hjust = 1, vjust = 0.5, size = 11),
        strip.text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "bottom",
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("")+ylab("")


hm.clean <- plot.heatmap.module.cor.contigs.lefse.neutral.sediment +
  theme(axis.title.y = element_blank(), axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none")

right <- ggplot(data = barchart.right, aes(x = module.from, y = n)) + 
  geom_bar(stat = "identity") + ylab("") + 
  scale_y_continuous(expand = c(0.01, 0), limits = c(0,95))+theme_classic()+
  geom_text(aes(label = n, y = n + 5), 
            angle = 0, hjust = 0.5, vjust = 0.5, size = 5, color = "red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y= element_blank(),
        legend.position="none")  +coord_flip()

up <- ggplot(data = barchart.up, aes(x = split, y = n)) + 
  geom_bar(stat = "identity") + ylab("") + 
  scale_y_continuous(expand = c(0.01, 0),  limits = c(0,68))+theme_classic()+ 
  geom_text(aes(label = n, y = n + 2), 
            angle = 360, hjust = 0.5, vjust = 0.5, size = 5, color = "red") +
  #geom_text(
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y= element_blank(),
        legend.position="none") 

cowplot::plot_grid(
  up, ggplot(), 
  hm.clean, right, 
  byrow = TRUE,
  nrow = 2, ncol = 2, rel_widths = c(1.5 ,0.3), rel_heights = c(0.3, 1.5)
)


##### c ######
ggplot(dat.figS7.heatmap.cor.subdivision.contigs, 
                                               aes(x = group, y = column, fill = mean)) +
  geom_tile(aes(x = group, y = column, ,
                fill = mean))+ 
  theme_classic()+ 
  ggh4x::facet_nested(~media + module, scales = "free", space = "free")+
  # scale_alpha_manual(values = c(0.5, 0.75, 1)) +
  scale_fill_gradient2(low = "grey100",  mid = "grey50", high = "red", midpoint = 0.2) +
  theme_classic()+
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      label <- gsub("above", "<span style='color:#A52A2A;'>above</span>", label)
      label <- gsub("below", "<span style='color:#29A6A6;'>below</span>", label)
      label <- gsub("neutral", "<span style='color:grey;'>neutral</span>", label)
      label
    })
  }) +
  geom_text(
    aes(x = group, y = column, label=letters), 
    color="black", size = 4) + 
  #facet_grid(~facet, scales = "free", space = "free")+
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 10),
        axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 11, angle = 45),
        axis.text.y = ggtext::element_markdown(hjust = 1, vjust = 0.5, size = 11),
        strip.text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "bottom",
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("")+ylab("")




