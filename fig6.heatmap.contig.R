library(dplyr)
library(stringr)
library(dunn.test)
library(multcompView)
### import arg ####
arg_gene_240927.tsv <- read.delim("G:/Bob/rproject2/Tibet_2023/Tibet2023/159_tibet_gene_240927/159_tibet_gene_240927/arg_gene_240927.tsv")
a <- arg_gene_240927.tsv
a <- a[,sapply(a, is.character)]
apply(a, 2, function(x) length(unique(x)))

a <- vf_gene_240927.tsv
a <- a[,sapply(a, is.character)]
apply(a, 2, function(x) length(unique(x)))

a <- mge_gene_240927.tsv
mge_gene_240927.tsv %>% group_by(Minor.mobileOG.Categories) %>% summarize(n = n()) %>% as.data.frame()
a <- a[,sapply(a, is.character)]
apply(a, 2, function(x) length(unique(x)))

#### #####
a <- arg_gene_240927.tsv
a <- a[,sapply(a, is.numeric)]
rownames(a) <- arg_gene_240927.tsv$orf1 
a$subtype <- arg_gene_240927.tsv$subtype
a <- aggregate(.~subtype, a, sum)
rownames(a) <- a$subtype

b <- data.frame(subtype = a$subtype)
b <- data.frame(name1 = str_split_fixed(string = rownames(a), pattern = "__", n = 2)[,1],
                name2 = str_split_fixed(string = rownames(a), pattern = "__", n = 2)[,2],
                subtype = rownames(a))
b$type <- arg_gene_240927.tsv$type[match(b$subtype, arg_gene_240927.tsv$subtype)]
b$alias <- paste("arg_gene", 1:746, sep = "_")
b
rownames(b) <- b$alias
rownames(a) <- b$alias
a$subtype <- NULL
b <- b[,c("type","subtype","alias")]
colnames(a) <- gsub(pattern = "RPKM", replacement = "", x = colnames(a))
phy.gene.arg <- phyloseq(otu_table(a, taxa_are_rows = TRUE), 
                         tax_table(as.matrix(b)), sample_data(meta.tab))


### import mge ####
mge_gene_240927.tsv <- read.delim("G:/Bob/rproject2/Tibet_2023/Tibet2023/159_tibet_gene_240927/159_tibet_gene_240927/mge_gene_240927.tsv")
a <- mge_gene_240927.tsv

a$Name[is.na(a$Name)] <- "NA"
a$Name2 <- ifelse(a$Name %in% c("NA:Keyword", "NA"), 
                  yes = paste(a$Name, a$mge_subtype, sep = "__"), 
                  no = a$Name)

mge_gene_240927.tsv <- a
a <- mge_gene_240927.tsv
a <- a[,sapply(a, is.numeric)]
a$orf1 <- mge_gene_240927.tsv$orf1
b <- a %>% group_by(orf1) %>% summarize(n = n()) %>% filter(n > 1) %>% as.data.frame

c <- mge_gene_240927.tsv[mge_gene_240927.tsv$orf1 %in%b$orf1,]
c <- c[!grepl(pattern = "Keyword", c$Name), ]

d <- mge_gene_240927.tsv[!mge_gene_240927.tsv$orf1 %in% b$orf1, ]
d <- rbind(c, d)
rownames(d) <- d$orf1
e <- d[,sapply(d, is.numeric)]
e$mge_subtype <- d$mge_subtype
e <- aggregate(.~mge_subtype, e, sum)
rownames(e) <- e$mge_subtype
e
f <- data.frame(mge_subtype = rownames(e))
rownames(f) <- f$mge_subtype

colnames(e) <- gsub(pattern = "RPKM", replacement = "", x = colnames(e))
rownames(f) <- f$mge_subtype
rownames(e) ==rownames(f)
e$mge_subtype <- NULL
f
phy.gene.mge <- phyloseq(otu_table(e, taxa_are_rows = TRUE), 
                         tax_table(as.matrix(f)), sample_data(meta.tab))


### import vf ####
vf_gene_240927.tsv <- read.delim("G:/Bob/rproject2/Tibet_2023/Tibet2023/159_tibet_gene_240927/159_tibet_gene_240927/vf_gene_240927.tsv")
vf_gene_240927.tsv$VFCID %>% unique
vf_gene_240927.tsv$VFCID %>% unique

a <- vf_gene_240927.tsv
a <- a[,sapply(a, is.character)]
apply(a, 2, function(x) length(unique(x)))
a$VF_short_name %>% unique

a <- vf_gene_240927.tsv
a <- a[,sapply(a, is.numeric)]
a$vf_geneID <- vf_gene_240927.tsv$vf_geneID 
a <- aggregate(.~vf_geneID, a, sum)
rownames(a) <- a$vf_geneID

f <- data.frame(vf_geneID = rownames(a))
f$VF_short_name <- vf_gene_240927.tsv$VF_short_name[match(f$vf_geneID, 
                                                          vf_gene_240927.tsv$vf_geneID)]
f$VF_name <- vf_gene_240927.tsv$VF_name[match(f$vf_geneID, 
                                              vf_gene_240927.tsv$vf_geneID)]
f$VFCID <- vf_gene_240927.tsv$VFCID[match(f$vf_geneID, 
                                          vf_gene_240927.tsv$vf_geneID)]
f$vf_Category <- vf_gene_240927.tsv$vf_Category[match(f$vf_geneID, 
                                                      vf_gene_240927.tsv$vf_geneID)]

colnames(a) <- gsub(pattern = "RPKM", replacement = "", x = colnames(a))
f <- f[,c("vf_Category", "VFCID", "VF_name","VF_short_name","vf_geneID")]
rownames(f) <- f$vf_geneID
rownames(a)
a$vf_geneID <- NULL
phy.gene.vf <- phyloseq(otu_table(a, taxa_are_rows = TRUE), 
                        tax_table(as.matrix(f)), sample_data(meta.tab))
phy.gene.vf
phy.gene.mge
phy.gene.arg

#### beta diversity #####
library(phyloseq)
a <- subset_samples(phy.gene.arg)
ord.a <- ordinate(physeq = a, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data1 <- data
data1$perc1 <- 34
data1$perc2 <- 14.7


a <- subset_samples(phy.gene.mge)
ord.a <- ordinate(physeq = a, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data2 <- data
data2$perc1 <- 26.9
data2$perc2 <- 12.8

ord.a <- ordinate(physeq = phy.gene.vf, method = "PCoA") ## nongut
plot_ordination(ordination = ord.a, physeq = a, 
                type = "samples", color="Tributary", shape = "Season")
data <- ord.a$vectors[,1:2] %>% as.data.frame
data$SampleID <- rownames(data)
data <- merge(data, meta.tab)
data3 <- data
data3$perc1 <- 47.1
data3$perc2 <- 12.5

data1$facet <- "ARG genes"
data2$facet <- "MGE genes"
data3$facet <- "VFG genes"
dim(data1);dim(data2);dim(data3)
data.pcoa.water.sd.gene <- rbind(data1, data2, data3)
data.pcoa.water.sd.gene$type <- factor(data.pcoa.water.sd.gene$type, levels = c("IF","WA","SD"), labels = c("Sewage", "Water", "Sediment"))

library(ggplot2)
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
fig.pcoa.gene

#### abundance heatmap #####
#### non sewage fig #####

a <- phy.gene.arg %>% otu_table %>% as.matrix %>% as.data.frame
a <- data.frame(arg = colSums(a))

b <- phy.gene.mge %>% otu_table %>% as.matrix %>% as.data.frame
b <- data.frame(mge = colSums(b))

c <- phy.gene.vf %>% otu_table %>% as.matrix %>% as.data.frame
c <- data.frame(vf = colSums(c))

d <- cbind(a, b, c) %>% as.data.frame()
d$SampleID <- rownames(d)
d <- merge(d, meta.tab[,c("SampleID", "type", "River3")])
colnames(d)[2:4] <- c("ARGs", "MGEs", "VFGs")
d$group <- paste(d$type, d$River3, sep = "_")
d$group %>% unique
d$group[grepl(pattern = "IF", x = d$group)] <- "Sewage"

cor.test(d[d$type =="WA", "ARGs"],
         d[d$type =="WA", "MGEs"])
cor.test(d[d$type =="WA", "ARGs"],
         d[d$type =="WA", "VFGs"])
cor.test(d[d$type =="SD", "ARGs"],
         d[d$type =="SD", "MGEs"])
cor.test(d[d$type =="SD", "ARGs"],
         d[d$type =="SD", "VFGs"])
plot(d[d$type =="IF", "ARGs"],
     d[d$type =="IF", "MGEs"])
plot(d[d$type =="IF", "ARGs"],
     d[d$type =="IF", "VFGs"])





#d <- d[,!colnames(d) %in% c("SampleID", "type", "River3")]
rownames(d) <- d$SampleID
e <- d[dat.SEM.water.new$SampleID, ]
cor.test(e$ARGs, dat.SEM.water.new.new$arg_total.abund)
plot(e$MGEs, dat.SEM.water.new.new$arg_total.abund)
plot(e$VFGs, dat.SEM.water.new.new$arg_total.abund)

library(ggpubr)
ggboxplot(data = d, x = "type", y = "ARGs", fill = "type") + 
  stat_compare_means(comparisons = list(c("IF","WA"), c("IF","SD"), c("WA","SD")))

d <- aggregate(.~group, d, mean) ### mean or sum
rownames(d) <- d$group
d$group <- NULL
e <- log10(d + 1)
e$group <- rownames(e)
e <- reshape::melt(e)

e$group <- gsub(pattern = "Parlung zangbo", replacement = "Parlung Zangbo", x = e$group)
e$group <- gsub(pattern = "Lhasa river", replacement = "Lhasa", x = e$group)

e$group %>% unique
e$group <- factor(e$group, levels = c("Sewage","WA_Chayu", "SD_Chayu", "WA_Lancang", "SD_Lancang", 
                                      "WA_Lhasa", "SD_Lhasa", "WA_Nianchu", "SD_Nianchu",
                                      "WA_Niyang", "SD_Niyang", "WA_Parlung Zangbo", "SD_Parlung Zangbo",
                                      "WA_Shiquan", "SD_Shiquan", "WA_Xiangquan", "SD_Xiangquan",
                                      "WA_Yarlung Zangbo", "SD_Yarlung Zangbo"))

e$River3 <- str_split_fixed(string = e$group, pattern = "_", n = 2)[,2]
e <- e[e$group != "Sewage", ]


#### contigs #####
#### non sewage test ####
a <- aggregate(.~pahogenicity, arg_mge_contig_abun_240902, sum)
a$pahogenicity <- paste("argmge", a$pahogenicity, sep = "_")
rownames(a) <- a$pahogenicity
a$pahogenicity <- NULL

b <- aggregate(.~pahogenicity, arg_vf_contig_abun_240902, sum)
b$pahogenicity <- paste("argvf", b$pahogenicity, sep = "_")
rownames(b) <- b$pahogenicity
b$pahogenicity <- NULL

c <- aggregate(.~pahogenicity, arg_mge_vf_contig_abun_240902, sum)
c$pahogenicity <- paste("argmgevf", c$pahogenicity, sep = "_")
rownames(c) <- c$pahogenicity
c$pahogenicity <- NULL
dim(a)
dim(b)
d <- cbind(t(a), t(b), t(c)) %>% as.data.frame()

d$SampleID <- rownames(d)
d$all_non_pathogen <- d$argmge_non_pathogen + d$argvf_non_pathogen + d$argmgevf_non_pathogen
d$all_pathogen <- d$argmge_pathogen + d$argvf_pathogen + d$argmgevf_pathogen

d <- merge(d, meta.tab[,c("SampleID", "type", "River3")])
d$group <- paste(d$type, d$River3, sep = "_")
d$group %>% unique
#d$group[grepl(pattern = "IF", x = d$group)] <- "Sewage"

d <- d[d$type != "IF", ]
f <- reshape::melt(d)
f <- f[,c("type", "River3","variable","value")]
f$group <- paste(f$River3, f$variable, sep = "__")
g <- f %>% group_by(group) %>% summarise(sum = sum(value)) %>% filter(sum > 0)
f <- f[f$group %in% g$group, colnames(f) %in% c("group","type","value")]
colnames(f)[1] <- "media"
f$media[f$media == "WA"] <- "water"
f$media[f$media == "SD"] <- "sediment"

f <- split(f,
           f = list(f$group))
f$Lancang__all_pathogen
# Function to perform Wilcoxon test and return results with significance
perform_wilcox_test <- function(df, table_name) {
  # Calculate means and medians for each media
  mean_water <- mean(df$value[df$media == "water"], na.rm = TRUE)
  mean_sediment <- mean(df$value[df$media == "sediment"], na.rm = TRUE)
  median_water <- median(df$value[df$media == "water"], na.rm = TRUE)
  median_sediment <- median(df$value[df$media == "sediment"], na.rm = TRUE)
  
  # Perform Wilcoxon rank-sum test (Mann-Whitney U test)
  test_result <- wilcox.test(value ~ media, data = df)
  
  # Get p-value
  p_value <- test_result$p.value
  
  # Initialize p-value columns
  p0.1_water <- p0.1_sediment <- NA
  p0.05_water <- p0.05_sediment <- NA
  
  # Check significance and assign letters
  if (p_value <= 0.05) {
    # Different bold letters for water and sediment when significant at p ≤ 0.05
    p0.05_water <- "A"  # Bold letter for water media
    p0.05_sediment <- "B"  # Bold letter for sediment media
  } else if (p_value <= 0.1) {
    # Different non-bold letters for water and sediment when significant at 0.05 < p ≤ 0.1
    p0.1_water <- "a"  # Non-bold letter for water media
    p0.1_sediment <- "b"  # Non-bold letter for sediment media
  }
  
  # Create the result table for each media
  result <- data.frame(
    Table = table_name,
    Media = c("water", "sediment"),
    Mean = c(mean_water, mean_sediment),
    Median = c(median_water, median_sediment),
    pval = rep(p_value, 2),  # Same p-value for both media
    p0.1 = c(p0.1_water, p0.1_sediment),
    p0.05 = c(p0.05_water, p0.05_sediment)
  )
  
  return(result)
}

# Apply the function to each dataframe in the list, including a table name for identification
result_list <- lapply(seq_along(f), function(i) {
  perform_wilcox_test(f[[i]], paste("table", i))
})
names(result_list) <- names(f)
# Combine the results into a single table
final_result <- bind_rows(result_list, .id = "Table")

final_result$gene <- str_split_fixed(final_result$Table, pattern = "__", n = 2)[,2]
final_result$River3 <- str_split_fixed(final_result$Table, pattern = "__", n = 2)[,1]
colnames(final_result)[2] <- "media"
final_result$label <- paste(final_result$p0.1, final_result$p0.05)
final_result$label <- gsub(pattern = "NA", x = final_result$label,replacement = "" )
final_result$label <- gsub(pattern = " ", x = final_result$label,replacement = "" )
colnames(final_result)[4]<- "value"
final_result$River3 %>% unique


####  sewage test ####
a <- aggregate(.~pahogenicity, arg_mge_contig_abun_240902, sum)
a$pahogenicity <- paste("argmge", a$pahogenicity, sep = "_")
rownames(a) <- a$pahogenicity
a$pahogenicity <- NULL

b <- aggregate(.~pahogenicity, arg_vf_contig_abun_240902, sum)
b$pahogenicity <- paste("argvf", b$pahogenicity, sep = "_")
rownames(b) <- b$pahogenicity
b$pahogenicity <- NULL

c <- aggregate(.~pahogenicity, arg_mge_vf_contig_abun_240902, sum)
c$pahogenicity <- paste("argmgevf", c$pahogenicity, sep = "_")
rownames(c) <- c$pahogenicity
c$pahogenicity <- NULL

d <- cbind(t(a), t(b), t(c)) %>% as.data.frame()

d$SampleID <- rownames(d)
d$all_non_pathogen <- d$argmge_non_pathogen + d$argvf_non_pathogen + d$argmgevf_non_pathogen
d$all_pathogen <- d$argmge_pathogen + d$argvf_pathogen + d$argmgevf_pathogen
d
d <- merge(d, meta.tab[,c("SampleID", "type", "River3")])

f <- reshape::melt(d)
colnames(f)[2] <- "media"
f$media[f$media == "WA"] <- "water"
f$media[f$media == "SD"] <- "sediment"
f$media[f$media == "IF"] <- "sewage"
f[1:5,]
f
colnames(f)[3] <- "group"

filtered_list  <- split(f, f = list(f$variable))
filtered_list$argmge_non_pathogen[1:5,]
results_with_letters <- lapply(names(filtered_list), function(name) {
  df <- filtered_list[[name]]
  
  # Kruskal-Wallis test
  kw_result <- kruskal.test(value ~ media, data = df)
  
  # Dunn's test
  dunn_result <- dunn.test(df$value, df$media,  ,kw = TRUE)
  
  # Extract pairwise p-values
  p_values <- dunn_result$P
  comparisons <- dunn_result$comparisons
  
  # Ensure unique media levels
  media_levels <- sort(unique(df$media))
  
  # Initialize p-value matrix
  p_matrix <- matrix(1, nrow = length(media_levels), ncol = length(media_levels),
                     dimnames = list(media_levels, media_levels))
  
  # Populate the matrix with pairwise p-values
  for (j in seq_along(comparisons)) {
    medias <- unlist(strsplit(comparisons[j], " - "))
    if (all(medias %in% media_levels)) {
      p_matrix[medias[1], medias[2]] <- p_values[j]
      p_matrix[medias[2], medias[1]] <- p_values[j]
    }
  }
  
  # Use multcompLetters to assign letters
  letters <- multcompLetters(p_matrix, threshold = 0.1)$Letters
  
  # Combine results into a dataframe
  letters_df <- data.frame(
    media = names(letters),
    letters = unname(letters),
    table_name = name, # Add the name of the dataframe
    stringsAsFactors = FALSE
  )
  
  # Return results
  list(
    kw_result = kw_result,
    dunn_result = dunn_result,
    letters = letters_df
  )
})

# Combine letter assignments across all datasets
letters_combined <- bind_rows(lapply(results_with_letters, function(result) result$letters))
dunn.test(filtered_list$all_pathogen$value, filtered_list$all_pathogen$media, kw = TRUE)
filtered_list$all_pathogen %>% group_by(media) %>% summarise(mean = mean(value), median = median(value))
a <- filtered_list$all_pathogen
ggboxplot(data = a[a$media!= "sewage", ], x = "media", y = "value") + stat_compare_means()

results_with_letters <- lapply(names(filtered_list), function(name) {
  df <- filtered_list[[name]]
  
  # Kruskal-Wallis test
  kw_result <- kruskal.test(value ~ media, data = df)
  
  # Dunn's test
  dunn_result <- dunn.test(df$value, df$media,  ,kw = TRUE)
  
  # Extract pairwise p-values
  p_values <- dunn_result$P
  comparisons <- dunn_result$comparisons
  
  # Ensure unique media levels
  media_levels <- sort(unique(df$media))
  
  # Initialize p-value matrix
  p_matrix <- matrix(1, nrow = length(media_levels), ncol = length(media_levels),
                     dimnames = list(media_levels, media_levels))
  
  # Populate the matrix with pairwise p-values
  for (j in seq_along(comparisons)) {
    medias <- unlist(strsplit(comparisons[j], " - "))
    if (all(medias %in% media_levels)) {
      p_matrix[medias[1], medias[2]] <- p_values[j]
      p_matrix[medias[2], medias[1]] <- p_values[j]
    }
  }
  
  # Use multcompLetters to assign letters
  letters <- multcompLetters(p_matrix, threshold = 0.05)$Letters
  
  # Combine results into a dataframe
  letters_df <- data.frame(
    media = names(letters),
    letters = unname(letters),
    table_name = name, # Add the name of the dataframe
    stringsAsFactors = FALSE
  )
  
  # Return results
  list(
    kw_result = kw_result,
    dunn_result = dunn_result,
    letters = letters_df
  )
})

# Combine letter assignments across all datasets
letters_combined2 <- bind_rows(lapply(results_with_letters, function(result) result$letters))

f <- reshape::melt(d)
f[1:5,]
colnames(f)[2] <- "media"
f$media[f$media == "WA"] <- "water"
f$media[f$media == "SD"] <- "sediment"
f$media[f$media == "IF"] <- "sewage"
f$River3<- "Whole"

f <- f %>% group_by(River3, media, variable) %>% 
  summarise(mean = mean(value), median = median(value)) %>% as.data.frame
colnames(f)[3] <- "table_name"
letters_combined[1:4,]

d <- merge(f, letters_combined, by = c("media", "table_name"))
colnames(d)[6] <- "p0.1"
d <- merge(d, letters_combined2, by = c("media","table_name"))
colnames(d)[7] <- "p0.05"
d$label <- d$p0.05
colnames(d)[4] <- "value"
final_result[1:2,]

colnames(final_result)[8] <- "table_name"
colnames(final_result)[4] <- "median"
colnames(final_result)[3] <- "value"
final_result[final_result$River3 == "Xiangquan", ]
final_result$p0.1[final_result$p0.1 == "a"] <- "c"
final_result$p0.1[final_result$p0.1 == "b"] <- "d"
final_result$p0.05 <- tolower(final_result$p0.05)
final_result$label <- ifelse(!is.na(final_result$p0.05), final_result$p0.05, final_result$p0.1)
d$table_name %>% unique
d %>% filter(table_name == "all_pathogen")
e <- rbind(d[,c("media","table_name","value","River3","p0.1","p0.05","label")],
           final_result[,c("media","table_name","value","River3","p0.1","p0.05","label")])

colnames(e)[2] <- "variable"
e$table_name %>% unique
e$variable <- factor(e$variable, labels = c("ARG-MGE non-pathogen contigs",
                                                                                        "ARG-MGE pathogen contigs",
                                                                                        "ARG-VFGs non-pathogen contigs",
                                                                                        "ARG-VFGs pathogen contigs",
                                                                                        "ARG-MGE-VFGs non-pathogen contigs",
                                                                                        "ARG-MGE-VFGs pathogen contigs",
                                                                  "Total non-pathogen contigs",
                                                                  "Total pathogen contigs"))

e$media <- factor(e$media, levels = c("sewage", "water", "sediment"))
e$River3 %>% unique
e$River3 <- factor(e$River3, levels = c("Whole","Xiangquan", "Shiquan", "Yarlung Zangbo", "Nianchu", 
                                        "Lhasa river", "Niyang","Parlung zangbo", "Lancang", "Chayu"),
                   labels = c("Whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))
library(ggplot2)
write.csv(e, "fig.data.heatmap.contig.abundance.csv")

e <- fig6.heatmap.contigs.pathogen.non.media$data
e$label <- tolower(e$label)
e$River3 <- factor(e$River3, 
                   levels = c("Whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))

chayu <- e[e$River3 == "Xiangquan", ]
chayu$River3 <- "Chayu"
chayu$value <- NA
chayu$label <- NA
e <- rbind(e, chayu)
e[e$River3 == "Chayu", ]
e$value[e$value == 0] <- NA
log10(e$value+1) %>% mean
log10(e$value+1) %>% min
log10(e$value+1) %>% max
f <- e[e$River3 == "Xiangquan", ]
f$River3 <- "Chayu"
f$value <- 0
f$label <- NA
g <- rbind(e,f)
g$variable %>% unique
h <- g[g$variable == "ARG-MGE-VFGs pathogen contigs" & g$value == 0 & g$media == "water", ]
i <- h
i$media <- "sediment"
j <- rbind(h,i,h,i,h,i,h,i)
j$River3 <- rep(c("Yarlung\nZangbo", "Niyang","Parlung\nZangbo", "Lancang"), each = 2)
k <- rbind(h,i)
k$variable <- "ARG-MGE-VFGs non-pathogen contigs"
k$River3 <-  "Lancang"
k <- rbind(k,k)
k$variable[3:4] <- c("ARG-MGE pathogen contigs")

g <- rbind(g,k)
g <- rbind(g,j)


fig6.heatmap.contigs.pathogen.non.media <- ggplot(g,
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
