library(ggplot2)

dat.figS4.heatmp.cor.all
ggplot(dat.figS4.heatmp.cor.all[dat.figS4.heatmp.cor.all$var2 %in% c("Sewage", "Sediment", "Water"), ],
       aes(x = var1, y = var2, fill = cor)) +
  geom_tile(color = "white")+ 
  ggh4x::facet_grid2(rows = vars(facet), cols = vars(media), scales = "free", space = "free")+
  #facet_grid(~facet, scales = "free",space = "free" )+
  #facet_wrap(~facet2)+
  #facet_grid(~facet*media, scales = "free",  space = "free") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_classic()+
  #facet_grid(~facet2, scales = "free", space = "free") +
  #ggforce::facet_col(vars(facet2), scales = "free", space = "free") + 
  geom_text(data = dat.figS4.heatmp.cor.all[dat.figS4.heatmp.cor.all$var2 %in% c("Sewage", "Sediment", "Water"), ], 
            aes(x = var1, y = var2, label=sig), color="black", size=4) + 
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 12),
        axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 12, angle = 90),
        axis.text.y = ggtext::element_markdown(hjust = 1, vjust = 0.5, size = 12),
        strip.text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("")+ylab("") 

