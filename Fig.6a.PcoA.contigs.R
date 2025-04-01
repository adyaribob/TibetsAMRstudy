library(phyloseq)
library(ggplot2)
a <- subset_samples(phy.contig.argmge)
ord.a <- ordinate(physeq = a, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data1 <- data
data1$perc1 <- 20.2
data1$perc2 <- 11.7


a <- subset_samples(phy.contig.argvf)
ord.a <- ordinate(physeq = a, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data2 <- data
data2$perc1 <- 12.7
data2$perc2 <- 8.8

ord.a <- ordinate(physeq = phy.contig.argmgevf, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data3 <- data
data3$perc1 <- 41.5
data3$perc2 <- 14.9

data1$facet <- "ARG-MGE\ncarrying contigs"
data2$facet <- "ARG-VFG\ncarrying contigs"
data3$facet <- "ARG-MGE-VFG\ncarrying contigs"

data.pcoa.water.sd.contigs <- rbind(data1, data2, data3)
data.pcoa.water.sd.contigs$type <- factor(data.pcoa.water.sd.contigs$type, levels = c("IF","WA","SD"), labels = c("Sewage", "Water", "Sediment"))
write.csv(data.pcoa.water.sd.contigs,"data.pcoa.water.sd.contigs.fig6.csv", col.names = FALSE)

fig.pcoa.contigs <- ggplot(data = data.pcoa.water.sd.contigs) +
  geom_point(mapping = aes(x = Axis.1, y = Axis.2, color = type), 
             size  = 3.4, alpha = 0.8)+ facet_wrap(~domain, scales = "free") + 
  scale_color_manual(values = c("Sediment" = "#747070ff",  "Water" = "#31b2e6ff",
                                "Sewage" = "#c88b5cff" ))+
  theme_classic()+ 
  facet_grid(~facet , scales = "free") +
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 10),
        legend.position = "none",
        axis.text.y = ggtext::element_markdown( hjust = 1, 
                                                vjust = 0.7, size = 8),
        axis.text.x = ggtext::element_markdown( hjust = 1, 
                                                vjust = 0.7, size = 11),
        strip.text = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
