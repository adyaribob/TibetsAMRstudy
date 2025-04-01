library(cowplot)
library(dplyr)
library(stringr)

#### combine random forest results  barchart ####
#### water ####

b <- dat.res.random.forest.genus.contig.water
b <- b %>% filter(X.IncMSE.pval < 0.05)
colnames(b)[1:2] <- c("IncMSE","pval")
b$IncNodePurity <- NULL; b$IncNodePurity.pval <- NULL

b$subdivision <- str_split_fixed(string = b$Predictor, pattern = "_",n = 2)[,1]
b$genus1 <- str_split_fixed(string = b$Predictor, pattern = "_",n = 3)[,2]
b$genus2 <- str_split_fixed(string = b$Predictor, pattern = "_",n = 4)[,3]
b$subdivision <- ifelse(b$genus1  == "X", paste(b$subdivision, "X",sep = "_"), b$subdivision)
b$genus1 <- ifelse(b$genus1  == "X", b$genus2, b$genus1)
b$genus2[b$genus1 == b$genus2] <- ""
b$genus2[grepl(pattern = "module", b$genus2)] <- ""
b$genus1 <- ifelse(b$genus2  == "", b$genus1, paste(b$genus1, b$genus2, sep = "_"))
b$genus2 <- NULL
b$module <- sub(".*_(.*)", "\\1", b$Predictor)
b$module <- gsub(pattern = "module", replacement = "", x = b$module)
b$predictor <- paste(b$module, b$subdivision, b$genus1, sep = "_")

b$sig <- ifelse(b$pval > 0.05, "", 
                ifelse(b$pval <= 0.05 &
                         b$pval >= 0.01, "*", 
                       ifelse(b$pval < 0.01 &
                                b$pval >=0.001, "**", 
                              ifelse(b$pval < 0.001, "***","***"))))
c <- b

#### sediment ####
b <- dat.res.random.forest.contig.genus.sediment
b <- b %>% filter(X.IncMSE.pval < 0.05)
b$X <- NULL

colnames(b)[1:2] <- c("IncMSE","pval")
b$IncNodePurity <- NULL; b$IncNodePurity.pval <- NULL
library(stringr)
b$subdivision <- str_split_fixed(string = b$Predictor, pattern = "_",n = 2)[,1]
b$genus1 <- str_split_fixed(string = b$Predictor, pattern = "_",n = 3)[,2]
b$genus2 <- str_split_fixed(string = b$Predictor, pattern = "_",n = 4)[,3]
b$subdivision <- ifelse(b$genus1  == "X", paste(b$subdivision, "X",sep = "_"), b$subdivision)
b$genus1 <- ifelse(b$genus1  == "X", b$genus2, b$genus1)
b$genus2[b$genus1 == b$genus2] <- ""
b$genus2[grepl(pattern = "module", b$genus2)] <- ""
b$genus1 <- ifelse(b$genus2  == "", b$genus1, paste(b$genus1, b$genus2, sep = "_"))
b$genus2 <- NULL
b$module <- sub(".*_(.*)", "\\1", b$Predictor)
b$module <- gsub(pattern = "module", replacement = "", x = b$module)
b$predictor <- paste(b$module, b$subdivision, b$genus1, sep = "_")
b[grepl(pattern = "IF", x = b$PredictedVariable),]
b$sig <- ifelse(b$pval > 0.05, "", 
                ifelse(b$pval <= 0.05 &
                         b$pval >= 0.01, "*", 
                       ifelse(b$pval < 0.01 &
                                b$pval >=0.001, "**", 
                              ifelse(b$pval < 0.001, "***","***"))))
c$media <- "water"
b$media <- "sediment"
b <- rbind(c,b)

b$media <- factor(b$media, levels = c("water","sediment"))
c <- b[grepl(pattern = "IF", x = b$PredictedVariable),]
b <- heatmap.randomforest.contig.genus$data
b$module %>% unique
b$module <- factor(b$module, level = c("1","2","3","5","6","9","10"))
b <- arrange(b, module)
b$predictor <- factor(b$predictor, levels = unique(b$predictor))
dat.heatmap.randomforest.contig.genus <- b
heatmap.randomforest.contig.genus <- ggplot(data = b[grepl(pattern = "IF", x = b$PredictedVariable),], 
       aes(x = factor(PredictedVariable), y = predictor, fill = IncMSE)) + geom_tile() + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
  theme_classic()+ 
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      label <- gsub("above", "<span style='color:#A52A2A;'>above</span>", label)
      label <- gsub("below", "<span style='color:#29A6A6;'>below</span>", label)
      label <- gsub("neutral", "<span style='color:grey;'>neutral</span>", label)
      label
    })
  }) +
  ggh4x::facet_grid2(cols = vars(media), scales = "free", space = "free") +
  geom_text(aes(x = PredictedVariable, y = predictor, label=sig), color="black", size=4) + 
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    #strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 13, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    #strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 

hm.clean <- heatmap.randomforest.contig.genus +
  theme(axis.title.y = element_blank(), axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none")


#### barchart #####
a <- random.forest.R2.genus.contigs.water
b <- result.R2.randomforest.sediment
a$media <- "water"
b$media <- "sediment"
b <- rbind(a,b)
b <- b[grepl("IF", x = b$PredictedVariable), ]
b$X <- NULL
b$r2[b$r2 < 0 ] <- 0

assign_significance <- function(p) {
  ifelse(
    p < 0.001, "***",
    ifelse(
      p < 0.01, "**",
      ifelse(
        p < 0.05, "*",
        ""
      )
    )
  )
}

b$sig <- assign_significance(b$pval)
b$media <- factor(b$media, levels = c("water","sediment"))
dat.barchart.randomforest.contig.genus<- b
up <- ggplot(data = b, aes(x = PredictedVariable, y = r2)) + 
  geom_bar(stat = "identity", fill = "red") + ylab("") + 
  ggh4x::facet_grid2(cols = vars(media), scales = "free", space = "free") +
  scale_y_continuous(expand = c(0.01, 0),  limits = c(0,1))+theme_classic()+ 
  geom_text(aes(label = sig, y = r2 + 0.05), 
            angle = 360, hjust = 0.5, vjust = 0.5, size = 5, color = "black") +
  #geom_text(
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), #axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), axis.ticks.y= element_blank(),
        legend.position="none") 
up
library(ggplot2)
plot_grid( up,heatmap.randomforest.contig.genus, align = "v", nrow = 2, rel_heights = c(0.3, 1.5))
