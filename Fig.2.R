library(ggplot2)
library(dplyr)

a <- subset_samples(phy.gene.arg.subtype)
ord.a <- ordinate(physeq = a, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data1 <- data
data1$perc1 <- 37.8
data1$perc2 <-  12.1

a <- subset_samples(phy.gene.mge.mobilleOGID)
ord.a <- ordinate(physeq = a, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data2 <- data
data2$perc1 <- 25.8
data2$perc2 <- 12.3

ord.a <- ordinate(physeq = phy.gene.vf_geneID, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data3 <- data
data3$perc1 <- 38.6
data3$perc2 <- 15.6
data1
data1$facet <- "ARG genes"
data2$facet <- "MGE genes"
data3$facet <- "VFG genes"
dim(data1);dim(data2);dim(data3)
data.pcoa.water.sd.gene <- rbind(data1, data2, data3)
data.pcoa.water.sd.gene$type <- factor(data.pcoa.water.sd.gene$type, levels = c("IF","WA","SD"), labels = c("Sewage", "Water", "Sediment"))
dat.fig.2.genes.pcoa <- data.pcoa.water.sd.gene

dat.fig.2.genes.pcoa <- read.csv("F:/Non_uper/IUE_submitted_data/Tibet/Data/dat.fig2.genes.pcoa.csv")
(fig.2.pcoa.gene <- ggplot(data = dat.fig.2.genes.pcoa) +
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
        panel.border = element_rect(color = "black", fill = NA, size = 1)))

tab <- read.csv("F:/Non_uper/IUE_submitted_data/Tibet/Data/dat.fig2.abundance.csv")
d <- read.csv("F:/Non_uper/IUE_submitted_data/Tibet/Data/dat.fig2.abundance.test.csv")
d$p0.05 <- tolower(d$p0.05)
d$label <- d$p0.05
tab$River3 %>% unique
tab$River3 <- factor(tab$River3, 
                     levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                                "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))
d$River3 <- factor(d$River3, 
                   levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))

tab$media <- factor(tab$media, levels = c("sewage","water","sediment"))
d$p0.1[is.na(d$p0.1)] <- ""
d[d$p0.1 == "a", "p0.1"] <- "c"
d[d$p0.1 == "b", "p0.1"] <- "d"
d$label <- ifelse(is.na(d$label), yes = d$p0.1, no = d$label)

(fig2.boxplotabundance.media <- ggplot(tab, 
                                      aes(x = media, y = log10(value+1), fill = media)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  #scale_x_discrete(expand = c(0.15, 0.05)) +
  theme_classic()+ ggh4x::facet_grid2(cols = vars(River3), rows = vars(gene), scales = "free_x", space = "free_x")+
  labs(
    y = "Gene abundances (RPKM)") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "sewage" = "#c88b5cff", "none" = "grey80"))+
  xlab("")+
  geom_text(data = d, aes(x = media, y = log10(value + (1*value)), label = label), 
            size = 4.5, inherit.aes = FALSE) +
  theme(
    axis.text.y.left  = element_text(size = 13, colour = "black"),
    strip.text.x = element_text(size = 13, colour = "black"), 
    strip.text.y = element_text(size = 13, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 14 ,face = "plain", hjust = 1, vjust = 1, angle = 45))
)

tab <- read.csv("F:/Non_uper/IUE_submitted_data/Tibet/Data/dat.fig2.diversity.csv")
d <- read.csv("F:/Non_uper/IUE_submitted_data/Tibet/Data/dat.fig2.diversity.test.csv")
d$p0.05 <- tolower(d$p0.05)
d$label <- d$p0.05
d$p0.1[is.na(d$p0.1)] <- ""
d[d$p0.1 == "a", "p0.1"] <- "c"
d[d$p0.1 == "b", "p0.1"] <- "d"
tab$River3 <- factor(tab$River3, 
                     levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                                "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))
d$River3 %>% unique
d$River3 <- factor(d$River3, 
                   levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))
tab$media <- factor(tab$media, levels = c("sewage","water","sediment"))
d$label <- ifelse(is.na(d$label), yes = d$p0.1, no = d$label)

(fig2.boxplotdiversity.media <- ggplot(tab, 
                                      aes(x = media, y = log10(value+1), fill = media)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  #scale_x_discrete(expand = c(0.15, 0.05)) +
  theme_classic()+ ggh4x::facet_grid2(cols = vars(River3), rows = vars(gene), 
                                      scales = "free_x", space = "free_x")+
  labs(
    y = "Gene diversities") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "sewage" = "#c88b5cff", "none" = "grey80"))+
  xlab("")+
  geom_text(data = d, aes(x = media, y = log10(value + (1*value)), label = label), 
            size = 4.5, inherit.aes = FALSE) +
  theme(
    axis.text.y.left  = element_text(size = 13, colour = "black"),
    strip.text.x = element_text(size = 13, colour = "black"), 
    strip.text.y = element_text(size = 13, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 14 ,face = "plain", hjust = 1, vjust = 1, angle = 45))
)
