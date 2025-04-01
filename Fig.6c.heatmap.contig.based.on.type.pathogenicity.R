
write.csv(dat.heatmap.contigs.type.pathogenicity.media.fig6c, "dat.heatmap.contigs.type.pathogenicity.media.fig6c.csv", col.names = FALSE)

fig6.heatmap.contigs.pathogen.non.media <- ggplot(dat.heatmap.contigs.type.pathogenicity.media.fig6c,
                                                  aes(x = media, y = variable, fill = log10(value+1))) +
  geom_tile(color = "white")+ #+ facet_wrap(~order, scales = "free_y")+
  scale_fill_gradient2(low = "grey90",  mid = "grey50", high = "red", midpoint = 1, na.value = "white") +
  theme_classic()+
  geom_text(data = g, 
            aes(x = media, y = variable, label=label), color="black", size=4.5) + 
  facet_grid(~River3, scales = "free", space = "free")+
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 10),
        axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 13, angle = 45),
        axis.text.y = ggtext::element_markdown(hjust = 1, vjust = 0.5, size = 13),
        strip.text = element_text(size = 13),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("")+ylab("")
fig6.heatmap.contigs.pathogen.non.media