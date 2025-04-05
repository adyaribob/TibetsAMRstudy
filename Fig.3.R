library(ggplot2)

a <- phy.gene.arg.subtype %>% otu_table %>% as.matrix %>% as.data.frame
b <- phy.gene.mge.mobilleOGID %>% otu_table %>% as.matrix %>% as.data.frame
rownames(b) <- paste("mge", sep = "_", rownames(b))
c <- phy.gene.vf_geneID %>% otu_table %>% as.matrix %>% as.data.frame


a$gene <- rownames(a); a$lefse <- data.lefse.gene.arg$Group[match(a$gene, data.lefse.gene.arg$group2)]
b$gene <- rownames(b); b$lefse <- data.lefse.gene.mge$Group[match(b$gene, data.lefse.gene.mge$group2)]
c$gene <- rownames(c); c$lefse <- data.lefse.gene.vf$Group[match(c$gene, data.lefse.gene.vf$group2)]

a$lefse[is.na(a$lefse)] <- "Nogroup"
b$lefse[is.na(b$lefse)] <- "Nogroup"
c$lefse[is.na(c$lefse)] <- "Nogroup"

a$lefse[a$lefse == "IF"] <- "Sewage"
b$lefse[b$lefse == "IF"] <- "Sewage"
c$lefse[c$lefse == "IF"] <- "Sewage"

a$lefse[a$lefse == "SD"] <- "Sediment"
b$lefse[b$lefse == "SD"] <- "Sediment"
c$lefse[c$lefse == "SD"] <- "Sediment"

a$lefse[a$lefse == "WA"] <- "Water"
b$lefse[b$lefse == "WA"] <- "Water"
c$lefse[c$lefse == "WA"] <- "Water"

a$domain <- "ARGs"
b$domain <- "MGEs"
c$domain <- "VFGs"

a$gene <- NULL
b$gene <- NULL
c$gene <- NULL

d <- rbind(a,b,c)
d <- d %>% group_by(lefse, domain) %>% summarize(n = n()) %>% as.data.frame()
dat.fig3.lefse.arg.mge.vf.barchart <- d

dat.fig3.lefse.arg.mge.vf.barchart <- read.csv("dat.fig3.lefse.arg.mge.vf.barchart.csv")
(fig.3a.barchart.gene.lefse <- ggplot(dat.fig3.lefse.arg.mge.vf.barchart[dat.fig3.lefse.arg.mge.vf.barchart$lefse != "Nogroup", ], 
                                  aes(x = lefse, y = n, fill = lefse)) +
  geom_bar(stat="identity") +
  ylab("Number of genes") + xlab("") + 
  scale_fill_manual(values = c("Sediment" = "#747070ff",  "Water" = "#31b2e6ff",
                               "Sewage" = "#c88b5cff" ))+
  ggh4x::facet_grid2(~domain, scales = "free", independent = "y")+
  theme_classic()+
  theme(
    axis.text.y.left  = element_text(size = 8, colour = "black"),
    strip.text.x = element_text(size = 10.5, colour = "black"), 
    strip.text.y = element_text(size = 10.5, colour = "black"), 
    axis.title.y = element_text(size = 10),
    strip.background = element_rect(size = 0.5),
    strip.text.x.top  = element_text(size = 9),
    strip.text.y.right =   element_text(size = 9),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 10.5 ,
                                           face = "plain", hjust = 1, vjust = 1, angle = 45)))

#### neutral model #####
dat.fig3.neutral.model.gene
dat.fig3.neutral.model.gene <- read.csv("dat.fig3.neutral.model.gene.csv")

(fig.3.neutral.genes <- ggplot(dat.fig3.neutral.model.gene, aes(x=log10(p), y=freq ))+
  geom_point(aes(color = lefse),alpha=0.5, size=1)+
  geom_line(aes(x=log10(p), y=freq.pred),
            color="blue",linetype="solid",size=0.7)+
  geom_line(aes(x=log10(p), y=pred.lwr),
            color="blue",linetype="dashed",size=0.7)+
  geom_line(aes(x=log10(p), y=pred.upr),
            color="blue",linetype="dashed",size=0.7)+
  scale_color_manual(name = "lefse", values = c("SD" = "#747070ff",  "WA" = "#31b2e6ff",
                                                "IF" = "#c88b5cff" , "None" = "purple1"))+ 
  ggh4x::facet_grid2(rows = vars(domain), cols = vars(media))+
  xlab("Log10 (Mean Relative abundance)") + ylab("Occurance frequency")+theme_classic()+
  theme(
    strip.text.x = element_text(size = 10.5, colour = "black"), 
    strip.text.y = element_text(size = 10.5, colour = "black"), 
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    strip.background = element_rect(size = 0.5),
    strip.text.x.top  = element_text(size = 9),
    strip.text.y.right =   element_text(size = 9),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 8 ,
                                           face = "plain", hjust = 1, vjust = 0.4))
)

### barchart ####
b <- filter(dat.fig3.neutral.model.gene, domain %in% c("ARGs", "MGEs", "VFGs")) %>% 
  group_by(media, domain, lefse, upper_neutral_under) %>% dplyr::summarize(n = n()) %>% as.data.frame
b$lefse[b$lefse == "IF"] <- "Sewage"
b$lefse[b$lefse == "WA"] <- "Water"
b$lefse[b$lefse == "SD"] <- "Sediment"
b$lefse <- factor(b$lefse, levels = c("Sewage","Water","Sediment","None"))
b$domain
b$n <- as.numeric(b$n)

(fig.3.neutral.genes.barchart <- ggplot(b, aes(x = lefse, y = n, fill = upper_neutral_under)) +
  geom_bar(stat = "identity", position = position_stack()) +
  ylab("Number of ARGs/VFGs") + xlab("") + 
  scale_fill_manual(name = "lefse", values = c("above" = "#A52A2A",  "neutral" = "grey",
                                               "below" = "#29A6A6" ))+
  
  ggh4x::facet_grid2(domain ~ media, scales = "free_y", independent = "y")+
  theme_classic()+
  theme(
    axis.text.y.left  = element_text(size = 8, colour = "black"),
    strip.text.x = element_text(size = 10.5, colour = "black"), 
    strip.text.y = element_text(size = 10.5, colour = "black"), 
    axis.title.y = element_text(size = 10),
    strip.background = element_rect(size = 0.5),
    strip.text.x.top  = element_text(size = 9),
    strip.text.y.right =   element_text(size = 9),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 10.5 ,
                                           face = "plain", hjust = 1, vjust = 1, angle = 45))
)

##### boxplot of abundance of genes based on lefse in water and sediment #######
dat.fig3.boxplot.abundance <- read.csv("dat.fig3.boxplot.abundance.csv")
dat.fig3.boxplot.abundance.sig.test <- read.csv("dat.fig3.boxplot.abundance.sig.test.csv")

ggplot(dat.fig3.boxplot.abundance, 
       aes(x = group, y = log10(value+1), fill = group)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  #scale_x_discrete(expand = c(0.15, 0.05)) +
  theme_classic()+ ggh4x::facet_grid2(cols = vars(River3), rows = vars(facet_row), scales = "free_x", space = "free_x")+
  labs(
    y = "Gene abundances (RPKM)") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "sewage" = "#c88b5cff", "none" = "grey80"))+
  xlab("")+
  geom_text(data = dat.fig3.boxplot.abundance.sig.test, 
            aes(x = group, y = log10(value + (1*value)), label = p0.0.5), 
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

