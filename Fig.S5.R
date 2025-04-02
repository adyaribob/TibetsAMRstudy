library(phyloseq)
library(dplyr)
library(Hmisc)

a <- read.delim("D:/rproject2/Tibet_2023/Tibet2023/feast_ARGgene_water/sinkwater_sourceinfsediment_allcontig2_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Sediment", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Water", replacement = "", x = rownames(a))

sediment <- rowSums(a[,grepl(pattern = "SD", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "SD|Unknown", x = colnames(a))])

names(sediment) == names(sewage)
feast.contribution.water <- data.frame(sediment = sediment, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)

feast.contribution.water <- feast.contribution.water[rownames(dat.contig.all.non.redundant.wa),]
dim(dat.contig.all.non.redundant.wa)
dim(feast.contribution.water)
rownames(dat.contig.all.non.redundant.wa) == rownames(feast.contribution.water)

dat.contig.all.non.redundant.wa
d <- dat.contig.all.non.redundant.wa
a <- d
a[a > 0 ] <- 1
dim(d)
d <- d[,colSums(a) >= 6]
dim(d)
d$SampleID <- NULL

# Compute correlation and p-value
cor_pval <- apply(d, 2, function(x) {
  cor_result <- cor.test(x, feast.contribution.water$sewage)  # Correlate each column with the single parameter
  c(correlation = cor_result$estimate, p_value = cor_result$p.value)
})

# Convert to a data frame for easier handling
result <- as.data.frame(t(cor_pval))
colnames(result) <- c("Correlation", "P_value")
result$name <- rownames(result)
result$lefse <- dat.res.lefse.contig.remove.redundant$Group[match(result$name, dat.res.lefse.contig.remove.redundant$group)]
result$media <-"water"

##### sediment  ######
a <- read.delim("D:/rproject2/Tibet_2023/Tibet2023/feast_ARGgene_water/sinksediment_sourceinfwater_allcontig_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Water", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Sediment", replacement = "", x = rownames(a))

water <- rowSums(a[,grepl(pattern = "WA", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "WA|Unknown", x = colnames(a))])

names(water) == names(sewage)
feast.contribution.water <- data.frame(water = water, sewage = sewage, 
                                       unknown = unknown)
feast.contribution.water$SampleID <- rownames(feast.contribution.water)
feast.contribution.water <- feast.contribution.water[rownames(dat.contig.all.non.redundant.sd),]
dim(dat.contig.all.non.redundant.sd)
dim(feast.contribution.water)
rownames(dat.contig.all.non.redundant.sd) == rownames(feast.contribution.water)


d <- dat.contig.all.non.redundant.sd
a <- d
a[a > 0 ] <- 1
dim(d)
d <- d[,colSums(a) >= 6]
dim(d)
d$SampleID <- NULL

# Compute correlation and p-value
cor_pval <- apply(d, 2, function(x) {
  cor_result <- cor.test(x, feast.contribution.water$sewage)  # Correlate each column with the single parameter
  c(correlation = cor_result$estimate, p_value = cor_result$p.value)
})

# Convert to a data frame for easier handling
result2 <- as.data.frame(t(cor_pval))
colnames(result2) <- c("Correlation", "P_value")
result2$name <- rownames(result2)
result2$lefse <- dat.res.lefse.contig.remove.redundant$Group[match(result2$name, dat.res.lefse.contig.remove.redundant$group)]
result2$media <-"sediment"

result$neutral <- dat.neutral.water.contig.all$upper_neutral_under[match(result$name, 
                                                                         rownames(dat.neutral.water.contig.all))]
result2$neutral <- dat.neutral.sediment.contig.all$upper_neutral_under[match(result2$name, 
                                                                             rownames(dat.neutral.sediment.contig.all))]

res <- rbind(result, result2)
res <- res %>% filter(!is.na(Correlation))
res$pathogen <- tax.contig$pathogen[match(res$name, tax.contig$contig)]
res[1:2,]
res %>% group_by(media, lefse, pathogen, neutral) %>% summarise(mean = mean(Correlation)) %>% as.data.frame
res %>% filter(P_value < 0.05) %>% group_by(media, lefse) %>% summarise(n = n(), mean = mean(Correlation)) %>% as.data.frame
res$lefse[is.na(res$lefse)] <- "None"
dat.cor.contig.feast.sewage.proportion <- res
res$lefse <- factor(res$lefse, levels = c("IF","WA","SD","None"), labels = c("Sewage","Water","Sediment","None"))
res$media <- factor(res$media, levels = c("water", "sediment"))

res$lefse %>% unique

library(ggpubr)
# Define the list of items
items <- c("Sewage", "Water", "Sediment", "None")

# Generate all non-redundant pairs
pairs <- combn(items, 2, simplify = FALSE)
plot.boxplot.cor.contig.lefse.sewage.feast$data
plot.boxplot.cor.contig.lefse.sewage.feast <- ggplot(res,  
                                                     aes(x = lefse, y = Correlation, fill = lefse)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  scale_x_discrete(expand = c(0.15, 0.05)) +
  theme_classic()+ ggh4x::facet_grid2(cols = vars(media))+
  #scale_y_continuous(
  #  limits = c(0, 1.4),          # Set y-axis range
  #  breaks = seq(0, 1, 0.2)    # Set labels at intervals of 0.2
  #)+   
  labs(
    y = "Pearson correlation coef.") + 
  scale_fill_manual(values = c("Sediment" = "#747070ff",  "Water" = "#31b2e6ff",
                               "Sewage" = "#c88b5cff", "None" = "grey80"))+
  xlab("")+ 
  stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
                     comparisons = list(c("Sewage", "Water"), c("Sewage", "Sediment"), 
                                        c("Sewage", "None"), c("Water", "Sediment"),
                                        c("Water","None","Sediment", "None")), 
                     symnum.args =list(
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),
                       symbols = c("***", "**", "*", ".","")
                     ))+
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))


a <- read.delim("D:/rproject2/Tibet_2023/Tibet2023/feast_ARGgene_water/sinkwater_sourceinfsediment_allcontig2_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Sediment", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Water", replacement = "", x = rownames(a))

sediment <- rowSums(a[,grepl(pattern = "SD", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "SD|Unknown", x = colnames(a))])

names(water) == names(sewage)
a <- data.frame(sediment = sediment, sewage = sewage, 
                unknown = unknown)
a$SampleID <- rownames(a)

a$media <- "water"

b <- read.delim("D:/rproject2/Tibet_2023/Tibet2023/feast_ARGgene_water/sinksediment_sourceinfwater_allcontig_source_contributions_matrix.txt", row.names=1)
colnames(b) <- gsub(pattern = "_Inf|_Water", replacement = "", x = colnames(b))
rownames(b) <- gsub(pattern = "_Sediment", replacement = "", x = rownames(b))

water <- rowSums(b[,grepl(pattern = "WA", x = colnames(b))])
unknown <- b[,grepl(pattern = "Unknown", x = colnames(b))]
sewage <- rowSums(b[,!grepl(pattern = "WA|Unknown", x = colnames(b))])

names(water) == names(sewage)
b <- data.frame(water = water, sewage = sewage, 
                unknown = unknown)
b$SampleID <- rownames(b)

a$media <- "water"
b$media <- "sediment"

items <- c("sewage", "water", "unknown")

# Generate all non-redundant pairs
pairs <- combn(items, 2, simplify = FALSE)
colnames(a)[1] <- "water"
b
c <- rbind(reshape::melt(a),
           reshape::melt(b))
c$media <- factor(c$media, levels = c("water","sediment"))
c$variable %>% unique
c$variable <- factor(c$variable, levels = c("sewage","water", "unknown"))
library(ggpubr)

plot.boxplot.feast.contig <- ggplot(c, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  scale_x_discrete(expand = c(0.15, 0.2)) +
  theme_classic()+ ggh4x::facet_grid2(cols = vars(media), scales = "free", independent = "y")+
  stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
                     comparisons = pairs ) + 
  labs(
    y = "FEAST source proportion") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "sewage" = "#c88b5cff", "unknown" = "grey80"))+
  xlab("")+ 
  theme_classic()+
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))

plot.boxplot.feast.contig2 <- ggplot(c[c$variable == "sewage", ], aes(x = media, y = value, fill = variable)) + 
  geom_boxplot(outlier.size = 0.5, linewidth = 0.1) + 
  scale_x_discrete(expand = c(0.3, 0.3)) +
  #theme_classic()+ ggh4x::facet_grid2(cols = vars(media), scales = "free", independent = "y")+
  stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4 ) + 
  labs(
    y = "FEAST source proportion") + 
  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
                               "sewage" = "#c88b5cff", "unknown" = "grey80"))+
  xlab("")+ 
  theme_classic()+
  theme(
    axis.text.y.left  = element_text(size = 12, colour = "black"),
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))

library(cowplot)
plot.boxplot.feast.contig.combined <- plot_grid(
  plot.boxplot.feast.contig, plot.boxplot.feast.contig2, 
  ncol = 2,            # Arrange in a single row
  rel_widths = c(2, 1) # Set the width ratio
)

c <- rbind(a[,c("sewage", "media", "total.contig")],
           b[,c("sewage", "media", "total.contig")])
c[1:5,]
c$media <- factor(c$media, levels = c("water","sediment"))
plot.cor.sewage.feast.total.contig <- ggplot(c, aes(x = sewage, y = total.contig)) + 
  geom_point(size = 2, alpha = 0.7, color = "#c88b5cff") +
  theme_classic()+ 
  ggh4x::facet_grid2(cols = vars(media),  scales = "free", independent = "y")+
  
  geom_smooth( method = "lm", se = TRUE) +  stat_cor()+
  
  #scale_y_continuous(
  #  limits = c(0, 1.4),          # Set y-axis range
  #  breaks = seq(0, 1, 0.2)    # Set labels at intervals of 0.2
  #)+   
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


plot_grid(
  plot.boxplot.feast.contig.combined, plot.boxplot.cor.contig.lefse.sewage.feast,
  ncol = 1, nrow = 2, rel_heights = c(0.5, 0.5)
)


a <- tax_table(subset_taxa(phy.contig.argmge, !contig %in% taxa_names(phy.contig.argmgevf))) %>% as.matrix %>% as.data.frame %>% filter(kindom != "") %>% mutate(contiggroup = "argmge")
b <- tax_table(subset_taxa(phy.contig.argvf, !contig %in% taxa_names(phy.contig.argmgevf))) %>% as.matrix %>% as.data.frame %>% filter(kindom != "") %>% mutate(contiggroup = "argvg")
c <- tax_table(phy.contig.argmgevf) %>% as.matrix %>% as.data.frame %>% filter(kindom != "")  %>% mutate(contiggroup = "argmgevf")
arg_mge_vf_contig_anno_240828$pathogen <- c$pathogen[match(arg_mge_vf_contig_anno_240828$contig, c$contig)]
arg_mge_vf_contig_anno_240828$species <- c$species[match(arg_mge_vf_contig_anno_240828$contig, c$contig)]
arg_mge_contig_anno_240828$pathogen <- a$pathogen[match(arg_mge_contig_anno_240828$contig, a$contig)]
arg_mge_contig_anno_240828$species <- a$species[match(arg_mge_contig_anno_240828$contig, a$contig)]
arg_vf_contig_anno_240828$pathogen <- b$pathogen[match(arg_vf_contig_anno_240828$contig, b$contig)]
arg_vf_contig_anno_240828$species <- b$species[match(arg_vf_contig_anno_240828$contig, b$contig)]

#####  #####
##### sig cor ####
#### water #####
a <- read.delim("D:/rproject2/Tibet_2023/Tibet2023/feast_ARGgene_water/sinkwater_sourceinfsediment_allcontig2_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Sediment", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Water", replacement = "", x = rownames(a))

sediment <- rowSums(a[,grepl(pattern = "SD", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "SD|Unknown", x = colnames(a))])

names(sediment) == names(sewage)
a <- data.frame(sediment = sediment, sewage = sewage, 
                unknown = unknown)
a$SampleID <- rownames(a)

a <- a[rownames(dat.contig.all.non.redundant.wa),]

rownames(dat.contig.all.non.redundant.wa) == rownames(a)
d <- dat.contig.all.non.redundant.wa

b <- d
b[b > 0 ] <- 1
dim(d)
d <- d[,colSums(b) >= 6]
dim(d)
d$SampleID <- NULL
rownames(d) == rownames(a)
table1 <- a
table2 <- d
table1$SampleID <- NULL

rownames(table1) == rownames(dat.contig.all.non.redundant.wa)
table1$total.contig <- rowSums(dat.contig.all.non.redundant.wa)

cor_matrix <- cor(as.matrix(table1), as.matrix(table2))
# Compute p-values
n <- nrow(table1)  # Number of observations
t_stats <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
p_matrix <- 2 * pt(-abs(t_stats), df = n - 2)

# Convert results to data frame
library(tidyr)
result <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Var1 = Var1, Var2 = Var2, correlation = Freq) %>%
  mutate(p_value = as.vector(p_matrix))

# View results
head(result)
result$neutral <- dat.neutral.water.contig.all$upper_neutral_under[match(result$Var2,
                                                                         rownames(dat.neutral.water.contig.all))]
result$lefse <- dat.res.lefse.contig.remove.redundant$Group[match(result$Var2, 
                                                                  dat.res.lefse.contig.remove.redundant$group)]
result$pathogen <- tax.contig$pathogen[match(result$Var2, tax.contig$contig)]

wide_df <- result %>%
  pivot_wider(
    names_from = Var1,
    values_from = c(correlation, p_value)
  )
wide_df$media <- "water"
list.cor.feast.with.totalabundance.and.feast  <- list(water = wide_df)


#### ####
#### sediment #####
a <- read.delim("D:/rproject2/Tibet_2023/Tibet2023/feast_ARGgene_water/sinksediment_sourceinfwater_allcontig_source_contributions_matrix.txt", row.names=1)
colnames(a) <- gsub(pattern = "_Inf|_Water", replacement = "", x = colnames(a))
rownames(a) <- gsub(pattern = "_Sediment", replacement = "", x = rownames(a))

water <- rowSums(a[,grepl(pattern = "WA", x = colnames(a))])
unknown <- a[,grepl(pattern = "Unknown", x = colnames(a))]
sewage <- rowSums(a[,!grepl(pattern = "WA|Unknown", x = colnames(a))])

names(water) == names(sewage)
a <- data.frame(water = water, sewage = sewage, 
                unknown = unknown)
a$SampleID <- rownames(a)
a <- a[rownames(a),]

a <- a[rownames(dat.contig.all.non.redundant.sd),]

rownames(dat.contig.all.non.redundant.sd) == rownames(a)
d <- dat.contig.all.non.redundant.sd
b <- d
b[b > 0 ] <- 1
dim(d)
d <- d[,colSums(b) >= 6]
dim(d)
rownames(d) == rownames(a)
table1 <- a
table2 <- d
table1$SampleID <- NULL
rownames(table1) == rownames(dat.contig.all.non.redundant.sd)
table1$total.contig <- rowSums(dat.contig.all.non.redundant.sd)

cor_matrix <- cor(as.matrix(table1), as.matrix(table2))
# Compute p-values
n <- nrow(table1)  # Number of observations
t_stats <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
p_matrix <- 2 * pt(-abs(t_stats), df = n - 2)

# Convert results to data frame
library(tidyr)
result <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Var1 = Var1, Var2 = Var2, correlation = Freq) %>%
  mutate(p_value = as.vector(p_matrix))

# View results
head(result)
result$neutral <- dat.neutral.sediment.contig.all$upper_neutral_under[match(result$Var2,
                                                                            rownames(dat.neutral.sediment.contig.all))]
result$lefse <- dat.res.lefse.contig.remove.redundant$Group[match(result$Var2, 
                                                                  dat.res.lefse.contig.remove.redundant$group)]
result$pathogen <- tax.contig$pathogen[match(result$Var2, tax.contig$contig)]
result[1:5,]
wide_df <- result %>%
  pivot_wider(
    names_from = Var1,
    values_from = c(correlation, p_value)
  )
wide_df$media <- "sediment"
list.cor.feast.with.totalabundance.and.feast$sediment <- wide_df
list.cor.feast.with.totalabundance.and.feast$sediment
saveRDS(list.cor.feast.with.totalabundance.and.feast, "list.cor.feast.with.totalabundance.and.feast.rds")
list.cor.feast.with.totalabundance.and.feast$water
