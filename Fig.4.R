

write.csv(dat.fig.module.tax.bac.euk.20240820, "dat.fig4.module.tax.bac.euk.20240820.csv", col.names = FALSE)
fig.modules.tax.composition <- ggplot(dat.fig.module.tax.bac.euk.20240820[dat.fig.module.tax.bac.euk.20240820$lefse2 != "Nogroup", ], 
                                      aes(x = lefse2, y = sum, fill = taxgroup)) + 
  geom_bar(stat = "identity",position = "stack") + 
  scale_fill_manual(values = dominant.phylum.subdivsion$color)+
  ggh4x::facet_grid2(vars(group2), vars(wtc.membership), 
                     scales = "free_y", independent = "y" )+
  labs(x = "", y = "Number of ASVs", fill = "Bacteria Phylum/\nBacteria Class/\nMicroeukaryotes subdivision")+
  theme_classic()+ 
  theme(
    axis.text.y.left  = element_text(size = 8, colour = "black"),
    strip.text.x = element_text(size = 10.5, colour = "black"), 
    strip.text.y = element_text(size = 9, colour = "black"), 
    axis.title.y = element_text(size = 10),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 10.5 ,face = "plain", hjust = 1, vjust = 0.6, angle = 90))

write.csv(dat.res.module.lefse, "dat.fig4.module.biomarker.strength.csv", col.names = FALSE)
fig.barchart.module.biomarker.strength <- ggplot(data = dat.res.module.lefse, aes(x = lefse, y = log10(score+1))) + 
  geom_bar(stat = "identity", aes(fill = lefse)) +
  #ggh4x::facet_nested(~media + module) + 
  theme_classic()+ 
  ggh4x::facet_grid2(vars(media), vars(module), scales = "free", independent = "y", space = "free")+
  scale_fill_manual(values = c("Sediment" = "#747070ff",  "Water" = "#31b2e6ff",
                               "Sewage" = "#c88b5cff" ))+
  labs(y = "Biomarker strength")+
  #geom_text(aes(x = group, y = Rsquare , label=sig), color="black", size=5) + 
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 7.5) +
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y.left  = element_text(size = 10, colour = "black"),
    strip.text.y = element_text(size = 9, colour = "black"), 
    axis.title.y = element_text(size = 10),
    legend.key = NULL,
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 9, angle = 45)  )


write.csv(dat.shared.network.sediment.water.with.inf.baceuk, "dat.fig4.shared.network.sediment.water.with.inf.baceuk.csv", col.names = FALSE)
fig.shared.edge.water.sediment.with.influent.barchart <- ggplot(data = dat.shared.network.sediment.water.with.inf.baceuk,
                                                                aes(x = module, y = rel)) + 
  geom_bar(stat = "identity", aes(fill = rel)) +
  scale_fill_gradient2(low = "white", high = "#c88b5cff", name = "Similarity to sewage")+
  ggh4x::facet_grid2(cols = vars(media), scales = "free", space = "free") + 
  theme_classic()+labs(x = "", y = "Number of shared edges\nwith sewage network (%)")+
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y.left  = element_text(size = 10, colour = "black"),
    strip.text.y = element_text(size = 9, colour = "black"), 
    axis.title.y = element_text(size = 10),
    legend.key = NULL,
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 10.5 ,face = "plain", hjust = 1, vjust = 1, angle = 45)
  )
