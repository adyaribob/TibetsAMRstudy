
fig.pcoa.gene <- ggplot(data = data.pcoa.water.sd.gene) +
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

tab <- read.csv("Fig1.abundance.csv")
d <- read.csv("Fig1.abundance.test.csv")
d$p0.05 <- tolower(d$p0.05)
d$label <- d$p0.05
tab$River3 %>% unique
tab$River3 <- factor(tab$River3, 
                     levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                                "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))
d$River3 %>% unique
d$River3 <- factor(d$River3, 
                   levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))

tab$media <- factor(tab$media, levels = c("sewage","water","sediment"))
d[d$p0.1 == "a", "p0.1"] <- "c"
d$p0.1 <- ifelse(d$p0.1 == "a", yes = "c", "d")
d$label <- ifelse(is.na(d$label), yes = d$p0.1, no = d$label)
fig1.boxplotabundance.media <- ggplot(tab, 
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


tab <- read.csv("Fig1.diversity.csv")
d <- read.csv("Fig1.diversity.test.csv")
d[d$label != "", ]
d$p0.05 <- tolower(d$p0.05)
d$label <- d$p0.05
tab$River3 %>% unique
tab$River3 <- factor(tab$River3, 
                     levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                                "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))
d$River3 %>% unique
d$River3 <- factor(d$River3, 
                   levels = c("whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))
tab$media <- factor(tab$media, levels = c("sewage","water","sediment"))
d$p0.1 <- ifelse(d$p0.1 == "a", yes = "c", "d")
d$label <- ifelse(is.na(d$label), yes = d$p0.1, no = d$label)

fig1.boxplotdiversity.media <- ggplot(tab, 
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
