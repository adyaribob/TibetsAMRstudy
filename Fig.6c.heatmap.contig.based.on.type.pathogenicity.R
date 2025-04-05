library(dplyr)
library(stringr)
library(dunn.test)
library(multcompView)
library(ggplot2)

### abundances of contigs based on their categories ####
arg_mge_contig_abun_240902
arg_vf_contig_abun_240902
arg_mge_vf_contig_abun_240902

#### contigs #####
#### water vs sediment: for each category of contigs based on genes they carry and pathogenicity of species  ####
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

d <- merge(d, meta.tab[,c("SampleID", "type", "River3")])
d$group <- paste(d$type, d$River3, sep = "_")
d$group %>% unique

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
    p0.05_water <- "a"  
    p0.05_sediment <- "b"  
  } else if (p_value <= 0.1) {
    p0.1_water <- "c"  
    p0.1_sediment <- "d"  
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
result

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

#### sewage vs. water vs. sediment: for each category of contigs based on genes they carry and pathogenicity of species  ####
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
d <- merge(d, meta.tab[,c("SampleID", "type", "River3")])

f <- reshape::melt(d)
colnames(f)[2] <- "media"
f$media[f$media == "WA"] <- "water"
f$media[f$media == "SD"] <- "sediment"
f$media[f$media == "IF"] <- "sewage"
colnames(f)[3] <- "group"

filtered_list  <- split(f, f = list(f$variable))

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
letters_combined <- bind_rows(lapply(results_with_letters, function(result) result$letters))

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
colnames(f)[2] <- "media"
f$media[f$media == "WA"] <- "water"
f$media[f$media == "SD"] <- "sediment"
f$media[f$media == "IF"] <- "sewage"
f$River3<- "Whole"

f <- f %>% group_by(River3, media, variable) %>% 
  summarise(mean = mean(value), median = median(value)) %>% as.data.frame
colnames(f)[3] <- "table_name"
letters_combined[1:4,]

g <- merge(f, letters_combined, by = c("media", "table_name"))
colnames(g)[6] <- "p0.1"
g <- merge(g, letters_combined2, by = c("media","table_name"))
colnames(g)[7] <- "p0.05"
g$label <- g$p0.05
colnames(g)[4] <- "value"

colnames(final_result)[8] <- "table_name"
colnames(final_result)[4] <- "median"
colnames(final_result)[3] <- "value"

final_result$label <- ifelse(!is.na(final_result$p0.05), final_result$p0.05, final_result$p0.1)
h <- rbind(g[,c("media","table_name","value","River3","p0.1","p0.05","label")],
           final_result[,c("media","table_name","value","River3","p0.1","p0.05","label")])

colnames(h)[2] <- "variable"
h$variable <- factor(h$variable, labels = c("ARG-MGE non-pathogen contigs",
                                            "ARG-MGE pathogen contigs",
                                            "ARG-VFGs non-pathogen contigs",
                                            "ARG-VFGs pathogen contigs",
                                            "ARG-MGE-VFGs non-pathogen contigs",
                                            "ARG-MGE-VFGs pathogen contigs",
                                            "Total non-pathogen contigs",
                                            "Total pathogen contigs"))

h$media <- factor(h$media, levels = c("sewage", "water", "sediment"))
h$River3 %>% unique
h$River3 <- factor(h$River3, levels = c("Whole","Xiangquan", "Shiquan", "Yarlung Zangbo", "Nianchu", 
                                        "Lhasa river", "Niyang","Parlung zangbo", "Lancang", "Chayu"),
                   labels = c("Whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))

chayu <- h[h$River3 == "Xiangquan", ]
chayu$River3 <- "Chayu"
chayu$value <- NA
chayu$label <- NA
chayu$p0.1 <- NA
h <- rbind(h, chayu)
h[h$River3 == "Chayu", ]
h$value[is.na(h$value)]<-0
dat.fig6c.heatmap.contigs.categories <- h
dat.fig6c.heatmap.contigs.categories <- read.csv("dat.fig6c.heatmap.contigs.categories.csv")

dat.fig6c.heatmap.contigs.categories$variable <- factor(dat.fig6c.heatmap.contigs.categories$variable, labels = c("ARG-MGE non-pathogen contigs",
                                            "ARG-MGE pathogen contigs",
                                            "ARG-VFGs non-pathogen contigs",
                                            "ARG-VFGs pathogen contigs",
                                            "ARG-MGE-VFGs non-pathogen contigs",
                                            "ARG-MGE-VFGs pathogen contigs",
                                            "Total non-pathogen contigs",
                                            "Total pathogen contigs"))

dat.fig6c.heatmap.contigs.categories$media <- factor(dat.fig6c.heatmap.contigs.categories$media, levels = c("sewage", "water", "sediment"))

dat.fig6c.heatmap.contigs.categories$River3 <- factor(dat.fig6c.heatmap.contigs.categories$River3, 
                                                     
                   levels = c("Whole", "Xiangquan", "Shiquan", "Yarlung\nZangbo", "Nianchu", 
                              "Lhasa", "Niyang","Parlung\nZangbo", "Lancang", "Chayu"))


(fig6c.heatmap.contigs.pathogen.non.media <- ggplot(dat.fig6c.heatmap.contigs.categories,
                                                  aes(x = media, y = variable, fill = log10(value+1))) +
  geom_tile(color = "white")+ #+ facet_wrap(~order, scales = "free_y")+
  scale_fill_gradient2(low = "grey90",  mid = "grey50", high = "red", midpoint = 1, na.value = "white") +
  theme_classic()+
  geom_text(data = dat.fig6c.heatmap.contigs.categories, 
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
  xlab("")+ylab(""))
