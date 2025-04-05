library(ggpubr)
""
##### water ####
##### MGE correlation to FEAST source track water ######
a <- read.delim("D:/.../sinkwater_sourceinfsediment_mge_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Sediment", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Water", replacement = "", x = rownames(a))

sediment <- rowSums(a[,grepl(pattern = "SD", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "SD|Unknown", x = colnames(a))])

names(sediment) == names(sewage)
feast.contribution.water <- data.frame(sediment = sediment, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)

result <- reshape::melt(feast.contribution.water)
result$gene <- "MGEs"

##### ARGs correlation to FEAST source track water ######
a <- read.delim("....sinkwater_sourceinfsediment_arg_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Sediment", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Water", replacement = "", x = rownames(a))

sediment <- rowSums(a[,grepl(pattern = "SD", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "SD|Unknown", x = colnames(a))])

names(sediment) == names(sewage)
feast.contribution.water <- data.frame(sediment = sediment, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)

result2 <- reshape::melt(feast.contribution.water)
result2$gene <- "ARGs"

##### VFGs correlation to FEAST source track water ######
a <- read.delim("....sinkwater_sourceinfsediment_vf_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Sediment", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Water", replacement = "", x = rownames(a))

sediment <- rowSums(a[,grepl(pattern = "SD", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "SD|Unknown", x = colnames(a))])

names(sediment) == names(sewage)
feast.contribution.water <- data.frame(sediment = sediment, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)

result3 <- reshape::melt(feast.contribution.water)
result3$gene <- "VFGs"

result4 <- rbind(result, result2, result3)
result4$media <- "water"

##### sediment ####
##### MGE correlation to FEAST source track water ######
a <- read.delim("....sinksediment_sourceinfwater_mge_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Water", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Sediment", replacement = "", x = rownames(a))

water <- rowSums(a[,grepl(pattern = "WA", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "WA|Unknown", x = colnames(a))])

names(water) == names(sewage)
feast.contribution.water <- data.frame(water = water, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)
result <- reshape::melt(feast.contribution.water)
result$gene <- "MGEs"

##### ARGs correlation to FEAST source track water ######
a <- read.delim("....sinksediment_sourceinfwater_arg_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Water", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Sediment", replacement = "", x = rownames(a))

water <- rowSums(a[,grepl(pattern = "WA", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "WA|Unknown", x = colnames(a))])

names(water) == names(sewage)
feast.contribution.water <- data.frame(water = water, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)
result2 <- reshape::melt(feast.contribution.water)
result2$gene <- "ARGs"

##### VFGs correlation to FEAST source track water ######
a <- read.delim("....sinksediment_sourceinfwater_vf_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Water", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Sediment", replacement = "", x = rownames(a))

water <- rowSums(a[,grepl(pattern = "WA", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "WA|Unknown", x = colnames(a))])

names(water) == names(sewage)
feast.contribution.water <- data.frame(water = water, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)

result3 <- reshape::melt(feast.contribution.water)
result3$gene <- "VFGs"

result5 <- rbind(result, result2, result3)
result5$media <- "sediment"

library(ggpubr)
result6 <- rbind(result4, result5)
result6[1:5,]
result6$variable <- factor(result6$variable, levels = c("sewage","water","sediment","unknown"), 
                           labels = c("Sewage","Water","Sediment","Unknown"))

## this 
result6$variable[result6$variable == "Sediment"] <- "Water"
# Define the list of items
items <- c("Sewage", "Water", "Unknown")

# Generate all non-redundant pairs
pairs <- combn(items, 2, simplify = FALSE)
result6$media <- factor(result6$media, levels = c("water", "sediment"))

dat.figS1.feast.source.proportion.boxplot.water <- read.csv("dat.figS1.feast.source.proportion.boxplot.water.csv")
ggplot(dat.figS1.feast.source.proportion.boxplot.water,  
       aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  scale_x_discrete(expand = c(0.15, 0.2)) +
  theme_classic()+ ggh4x::facet_grid2(rows = vars(media), cols = vars(gene), scales = "free", independent = "y")+
  stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
                     comparisons = pairs )+
  labs(
    y = "FEAST source proportion") + 
  scale_fill_manual(values = c("Sediment" = "#747070ff",  "Water" = "#31b2e6ff",
                               "Sewage" = "#c88b5cff", "Unknown" = "grey80"))+
  xlab("")+ 
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))+
  theme_classic()

ggplot(filter(dat.figS1.feast.source.proportion.boxplot.water, variable == "Sewage"),  
       aes(x = media, y = value), fill = "#c88b5cff") + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1,fill = "#c88b5cff") + 
  scale_x_discrete(expand = c(0.3, 0.3)) +
  theme_classic()+ ggh4x::facet_grid2(cols = vars(gene), scales = "free", independent = "y")+
  stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4) +
  labs(
    y = "FEAST source proportion of sewage", x = "")  +
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))

ggplot(dat.figS1.corplot.abundance.sewage.source.proportion[dat.figS1.corplot.abundance.sewage.source.proportion$variable == "Sewage", ], 
       aes(x = value, y = total, color = variable)) + 
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Sediment" = "#747070ff",  "Water" = "#31b2e6ff",
                                "Sewage" = "#c88b5cff", "unknown" = "grey80"))+
  theme_classic()+ 
  ggh4x::facet_grid2(rows = vars(media), cols = vars(gene), scales = "free", independent = "y")+
  
  geom_smooth(aes(color = variable), method = "lm", se = TRUE) +  stat_cor()+
  labs(
    y = "Total abundance (RPKM)") + 
  xlab("proportion of sewage source")+ 
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))

dat.figS1.boxplot.cor.sewagefeast.with.sewagesourceproportion <- read.csv("dat.figS1.boxplot.cor.sewagefeast.with.sewagesourceproportion.csv")

ggplot(dat.figS1.boxplot.cor.sewagefeast.with.sewagesourceproportion[dat.figS1.boxplot.cor.sewagefeast.with.sewagesourceproportion$from == "sewage.FEAST" , ],  
       aes(x = lefse, y = Correlation, fill = lefse)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  scale_x_discrete(expand = c(0.15, 0.05)) +
  theme_classic()+ ggh4x::facet_grid2(rows = vars(media), cols = vars(gene))+
  labs(
    y = "Pearson correlation coef.") + 
  scale_fill_manual(values = c("Sediment" = "#747070ff",  "Water" = "#31b2e6ff",
                               "Sewage" = "#c88b5cff", "None" = "grey80"))+
  xlab("")+ 
  stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
                     comparisons = pairs )+
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))
