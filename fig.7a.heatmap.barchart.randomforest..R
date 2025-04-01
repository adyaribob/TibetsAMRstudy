
dat.heatmap.randomforest.gene.contig.modules <- b
write.csv(dat.heatmap.randomforest.gene.contig.modules, "dat.heatmap.randomforest.gene.contig.modules.csv", col.names = FALSE)

fig.random.forest.heatmap <- ggplot(data = dat.heatmap.randomforest.gene.contig.modules, aes(x = factor(group), y = module, fill = IncMSE)) + geom_tile() + 
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
    axis.text.x = ggtext::element_markdown(size = 10.5 ,face = "plain", hjust = 1, vjust = 1, angle = 45))


##### random forest barchart #######

a <- dat.res.random.forest.water.contig.R2
a$group %>% unique
a <- a[a$group %in% c("contig_argvf", "contig_argmge", "contig_argmgevf"), ]
a$media <- "Water"
b <- dat.res.random.forest.sediment.contig.R2
b$group %>% unique
b <- b[b$group %in% c("contig_argvf", "contig_argmge", "contig_argmgevf"), ]
b$media <- "Sediment"


dat.res.random.forest.water.gene.arg.mge.vf.lefse.R2$media <- "Water"
dat.res.random.forest.sediment.gene.arg.mge.vf.lefse.R2$media <- "Sediment"

b <- rbind(dat.res.random.forest.water.gene.arg.mge.vf.lefse.R2,
           dat.res.random.forest.sediment.gene.arg.mge.vf.lefse.R2,
           a,
           b)

c$group %>% unique
c$group <- factor(c$group, levels = c("lefse_ARGs_IF", "lefse_VFGs_IF", "lefse_MGEs_IF",
                                      "lefse_ARGs_WA", "lefse_VFGs_WA", "lefse_MGEs_WA",
                                      "lefse_ARGs_SD", "lefse_VFGs_SD", "lefse_MGEs_SD",
                                      "contig_argmge", "contig_argvf", "contig_argmgevf"), 
                  labels = c("Sewage ARGs", "Sewage VFGs", "Sewage MGEs",
                             "Water ARGs", "Water VFGs", "Water MGEs",
                             "Sediment ARGs", "Sediment VFGs", "Sediment MGEs",
                             "ARG-MGE contigs", "ARG-VFG contigs", "ARG-MGE-VFG contigs"))

c$sig <- ifelse(c$pval > 0.05, "", 
                ifelse(c$pval <= 0.05 &
                         c$pval >= 0.01, "*", 
                       ifelse(c$pval < 0.01 &
                                c$pval >=0.001, "**", 
                              ifelse(c$pval < 0.001, "***","***"))))
c$media <- factor(c$media, levels = c("Water", "Sediment"))
c$Rsquare <- ifelse(c$Rsquare < 0, yes = 0, no = c$Rsquare)
dat.barchart.randomforest.gene.contig.modules<- c

write.csv(dat.barchart.randomforest.gene.contig.modules, "dat.barchart.randomforest.gene.contig.modules.csv", col.names = FALSE)

fig.random.forest.barchart <- ggplot(data = dat.barchart.randomforest.gene.contig.modules, aes(x = group, y = Rsquare)) + 
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

library(cowplot)
plot_grid(fig.random.forest.barchart, fig.random.forest.heatmap, 
          nrow = 2, cols = 1, align = "v", rel_widths = c(3,3), 
          rel_heights = c(1.2, 3))

