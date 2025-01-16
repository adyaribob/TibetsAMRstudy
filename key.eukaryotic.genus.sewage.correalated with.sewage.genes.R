library(phyloseq)
library(dplyr)
library(Hmisc)

### water #####
# 1. choose modules that are strongly associted with sewage
# 2. Sewage ARGS based on neutral partition (above, within, below) and pathogenicity (pathogen, non-pathogen)
# 3. average correlation between each sub-division of modules with the ARGs (according to above paritition)
# 4. for the n < 10 number of correlations between group of ARGs and subdivision of a modules, exculde them to get a rather reasonable n 
g4 <- dat.cor.euk.module.contig.argmge
g4$module.from <- dat.modularity.water.bac.euk.waktrap$wtc.membership[match(g4$from, 
                                                                            dat.modularity.water.bac.euk.waktrap$wtc.names)]
g4$pathogen <- tax.contig$pathogen[match(g4$to, tax.contig$contig)]
g4$lefse <- dat.res.lefse.contig.remove.redundant$Group[match(g4$to, dat.res.lefse.contig.remove.redundant$group)]
g4$neutral <- dat.neutral.water.contig.all$upper_neutral_under[match(g4$to, rownames(dat.neutral.water.contig.all))]
g4$subdivision <- tax$Subdivision[match(g4$from, tax$ASV)]
g4$genus <- tax$Genus[match(g4$from, tax$ASV)]
g4$family <- tax$Family[match(g4$from, tax$ASV)]
g4$order <- tax$Order[match(g4$from, tax$ASV)]
g4$tax <- paste(g4$subdivision, g4$genus, sep = ";")
colnames(g4)[3] <- "value"

g4$split <- paste(g4$subdivision, g4$genus,  g4$module.from, sep = "_")
g4$group <- paste(g4$neutral, g4$pathogen,  sep = "_")

g5 <- g4  %>% filter(lefse == "IF")
g5[1:5,]
## includes only module strongly related to influent 
check <- g5[g5$module.from %in% c(2,3,6,9,10), c("value", "group", "split")] %>% group_by(group, split) %>% 
  summarise(n = n()) %>%  as.data.frame %>% filter(n >= 10) 

### check that there is enough n to compare to others, 
# n = number correlation between ASVs within certain subdivision of certain module to group of ARGs (neutral,pathogenicity)

checklist <- split(check[,c("n", "group", "split")], f = list(check$split))
df_list <- split(g5[,c("value", "group", "split")], f = list(g5$split))

g5[g5$subdivision == "Gyrista" & g5$module.from == 9 & g5$neutral == "above" & g5$pathogen == "pathogen", "from"] %>% unique
g5[g5$subdivision == "Gyrista" & g5$module.from == 9 & g5$neutral == "above" & g5$pathogen == "pathogen", "to"] %>% unique

# Perform Dunn's test for each dataframe and assign letters
library(dunn.test)
library(multcompView)
library(dplyr)
library(stringr)

names(df_list)
names(checklist)

### check that there is no single group for each table
df_list <- df_list[sapply(df_list, function(df) length(unique(df[["group"]])) > 1)]
checklist <- checklist[sapply(checklist, function(df) length(unique(df[["group"]])) > 1)]

## only includes group in checklist (specific modules associated with sewage)
df_list <- df_list[names(checklist)]
df1 <- df_list
df2 <- checklist

# Filter based on unique values in corresponding data frames
filtered_list <- mapply(
  FUN = function(df1, df2) {
    df1[df1$group %in% unique(df2$group), ]
  },
  df1, df2,
  SIMPLIFY = FALSE
)
# View results

dunn_result <- dunn.test(filtered_list$IF_above_pathogen_8$value, 
                         filtered_list$IF_above_pathogen_8$group)
dunn_result

group_levels <- sort(unique(filtered_list$IF_above_pathogen_8$group))

results_with_letters <- lapply(names(filtered_list), function(name) {
  df <- filtered_list[[name]]
  
  # Kruskal-Wallis test
  kw_result <- kruskal.test(value ~ group, data = df)
  
  # Dunn's test
  dunn_result <- dunn.test(df$value, df$group,  ,kw = TRUE)
  
  # Extract pairwise p-values
  p_values <- dunn_result$P
  comparisons <- dunn_result$comparisons
  
  # Ensure unique group levels
  group_levels <- sort(unique(df$group))
  
  # Initialize p-value matrix
  p_matrix <- matrix(1, nrow = length(group_levels), ncol = length(group_levels),
                     dimnames = list(group_levels, group_levels))
  
  # Populate the matrix with pairwise p-values
  for (j in seq_along(comparisons)) {
    groups <- unlist(strsplit(comparisons[j], " - "))
    if (all(groups %in% group_levels)) {
      p_matrix[groups[1], groups[2]] <- p_values[j]
      p_matrix[groups[2], groups[1]] <- p_values[j]
    }
  }
  
  # Use multcompLetters to assign letters
  letters <- multcompLetters(p_matrix, threshold = 0.05)$Letters
  
  # Combine results into a dataframe
  letters_df <- data.frame(
    group = names(letters),
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
colnames(letters_combined)[3] <- "split"

SE <- function(x) {
  sd(x) / sqrt(length(x))
}

list.cor.feast.with.totalabundance.and.feast <- readRDS("list.cor.feast.with.totalabundance.and.feast.rds")

a <- list.cor.feast.with.totalabundance.and.feast$water
a[a$Var2 %in% unique(g5[g5$module.from == 5 & g5$subdivision == "Fungi" & g5$group == "neutral_pathogen", "to"]),
] %>% as.data.frame()
a[a$Var2 %in% unique(g5[g5$module.from == 2 , "to"]), ] %>% as.data.frame

### the group of genes and division with more than or equal to 10 number of correlation proceed, 
g6 <- g5 %>% dplyr::group_by(group, split) %>% 
  dplyr::summarise(n = n(), mean = mean(value), sd = sd(value), se = SE(value)) %>% as.data.frame()

### the dim of letters combined should be low because it does not have all modules, as in g6
## g6 has all modules but with m >= 10
dim(letters_combined);dim(g6);dim(merged_df)
merged_df <- merge(x = letters_combined, y = g6, by = c("group", "split"), all.x = TRUE)
merged_df <- merged_df[!is.na(merged_df$n), ]

merged_df$neutral <- g4$neutral[match(merged_df$split, g4$split)]
merged_df$pathogen <- g4$pathogen[match(merged_df$split, g4$split)]
merged_df$module <- g4$module.from[match(merged_df$split, g4$split)]
merged_df$rows <- paste(merged_df$neutral, merged_df$pathogen, sep = "\n")
merged_df$column <- str_split_fixed(merged_df$split, pattern = "_", n = 2)[,1]
merged_df$facet <- str_split_fixed(merged_df$split, pattern = "_", n = 2)[,2]
merged_df$rows <- gsub(pattern = "_", replacement = " ", x = merged_df$rows)
merged_df$group <- gsub(pattern = "_", replacement = " ", x = merged_df$group)
merged_df_water.genus <- merged_df

arrange(merged_df, desc(mean))

dim(merged_df_water)
dim(merged_df)
merged_df[1:5,]
merged_df_water[1:2,]

### sediment #####
# 1. choose all module as no sig. modules associated with sewage
# 2. Sewage ARGS based on neutral partition (above, within, below) and pathogenicity (pathogen, non-pathogen)
# 3. average correlation between each sub-division of modules with the ARGs (according to above paritition)
# 4. for the n < 10 number of correlations between group of ARGs and subdivision of a modules, exculde them to get a rather reasonable n 

g4 <- dat.cor.euk.module.contig.argmge.sediment

g4$module.from <- dat.modularity.sediment.bac.euk.waktrap$wtc.membership[match(g4$from, 
                                                                               dat.modularity.sediment.bac.euk.waktrap$wtc.names)]
g4$pathogen <- tax.contig$pathogen[match(g4$to, tax.contig$contig)]
g4$lefse <- dat.res.lefse.contig.remove.redundant$Group[match(g4$to, dat.res.lefse.contig.remove.redundant$group)]
g4$neutral <- dat.neutral.sediment.contig.all$upper_neutral_under[match(g4$to, rownames(dat.neutral.sediment.contig.all))]
g4$subdivision <- tax$Subdivision[match(g4$from, tax$ASV)]
g4$genus <- tax$Genus[match(g4$from, tax$ASV)]
g4$family <- tax$Family[match(g4$from, tax$ASV)]
g4$order <- tax$Order[match(g4$from, tax$ASV)]
g4$tax <- paste(g4$subdivision, g4$genus, sep = ";")

colnames(g4)[3] <- "value"

g4$split <- paste(g4$subdivision, g4$genus,  g4$module.from, sep = "_")
g4$group <- paste(g4$lefse, g4$neutral, g4$pathogen,  sep = "_")
g4[1:5,]
#g5 <- g4  %>% filter(genus != "" & lefse == "IF" & value > 0)
g5 <- g4  #%>% filter(value > 0)

check <- g5[,c("value", "group", "split")] %>% group_by(group, split) %>% 
  summarise(n = n()) %>%  as.data.frame 
#### select these three groups because others dont have too many contigs
check <- check[check$group %in% c("SD_above_non_pathogen", "SD_above_pathogen", "SD_neutral_non_pathogen"), ] %>% filter(n >= 10) 

checklist <- split(check[,c("n", "group", "split")], f = list(check$split))
df_list <- split(g5[,c("value", "group", "split")], f = list(g5$split))

# Perform Dunn's test for each dataframe and assign letters
library(dunn.test)
library(multcompView)
library(dplyr)
df_list <- df_list[sapply(df_list, function(df) length(unique(df[["group"]])) > 1)]
checklist <- checklist[sapply(checklist, function(df) length(unique(df[["group"]])) > 1)]

names(checklist)
names(df_list)

df_list <- df_list[names(checklist)]

df1 <- df_list
df2 <- checklist

#Filter based on unique values in corresponding data frames
filtered_list <- mapply(
  FUN = function(df1, df2) {
    df1[df1$group %in% unique(df2$group), ]
  },
  df1, df2,
  SIMPLIFY = FALSE
)
# View results


results_with_letters <- lapply(names(filtered_list), function(name) {
  df <- filtered_list[[name]]
  
  # Kruskal-Wallis test
  kw_result <- kruskal.test(value ~ group, data = df)
  
  # Dunn's test
  dunn_result <- dunn.test(df$value, df$group,  ,kw = TRUE)
  
  # Extract pairwise p-values
  p_values <- dunn_result$P
  comparisons <- dunn_result$comparisons
  
  # Ensure unique group levels
  group_levels <- sort(unique(df$group))
  
  # Initialize p-value matrix
  p_matrix <- matrix(1, nrow = length(group_levels), ncol = length(group_levels),
                     dimnames = list(group_levels, group_levels))
  
  # Populate the matrix with pairwise p-values
  for (j in seq_along(comparisons)) {
    groups <- unlist(strsplit(comparisons[j], " - "))
    if (all(groups %in% group_levels)) {
      p_matrix[groups[1], groups[2]] <- p_values[j]
      p_matrix[groups[2], groups[1]] <- p_values[j]
    }
  }
  
  # Use multcompLetters to assign letters
  letters <- multcompLetters(p_matrix, threshold = 0.05)$Letters
  
  # Combine results into a dataframe
  letters_df <- data.frame(
    group = names(letters),
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
colnames(letters_combined)[3] <- "split"

SE <- function(x) {
  sd(x) / sqrt(length(x))
}
g5[1:5,]
g6 <- g5 %>% dplyr::group_by(group, split) %>% dplyr::summarise(n = n(), mean = mean(value), 
                                                                sd = sd(value), se = SE(value)) %>% as.data.frame()
g6
merged_df <- merge(x = letters_combined, y = g6, by = c("group", "split"), all.x = TRUE)
dim(merged_df);dim(letters_combined)

merged_df <- merged_df[!is.na(merged_df$group), ]
merged_df$group %>% unique
merged_df$neutral <- g4$neutral[match(merged_df$split, g4$split)]
merged_df$pathogen <- g4$pathogen[match(merged_df$split, g4$split)]
merged_df$module <- g4$module.from[match(merged_df$split, g4$split)]
merged_df$rows <- paste(merged_df$neutral, merged_df$pathogen, sep = "\n")
merged_df$column <- str_split_fixed(merged_df$split, pattern = "_", n = 2)[,1]
merged_df$facet <- str_split_fixed(merged_df$split, pattern = "_", n = 2)[,2]
merged_df$rows <- gsub(pattern = "_", replacement = " ", x = merged_df$rows)
merged_df$group <- gsub(pattern = "_", replacement = " ", x = merged_df$group)
merged_df[1:5,]
merged_df_sediment.genus <- merged_df
merged_df_water$x <- paste0(merged_df_water$group, "(", merged_df_water$n, ")")
merged_df_sediment$x <- paste0(merged_df_sediment$group, "(", merged_df_sediment$n, ")")
merged_df_water$media <- NULL
merged_df_sediment$media <- NULL

merged_df_sediment$media <- "sediment"
merged_df_water$media <- "water"
a <- rbind(merged_df_water,merged_df_sediment)
a$media <- factor(a$media, levels = c("water","sediment"))

fig8.heatmap.cor.subdivision.contigs <- ggplot(a, 
                                               aes(x = group, y = column, fill = mean)) +
  geom_tile(aes(x = group, y = column, ,
                fill = mean))+ 
  theme_classic()+ 
  ggh4x::facet_nested(~media + module, scales = "free", space = "free")+
  # scale_alpha_manual(values = c(0.5, 0.75, 1)) +
  scale_fill_gradient2(low = "grey100",  mid = "grey50", high = "red", midpoint = 0.2) +
  theme_classic()+
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      label <- gsub("above", "<span style='color:#A52A2A;'>above</span>", label)
      label <- gsub("below", "<span style='color:#29A6A6;'>below</span>", label)
      label <- gsub("neutral", "<span style='color:grey;'>neutral</span>", label)
      label
    })
  }) +
  geom_text(
    aes(x = group, y = column, label=letters), 
    color="black", size = 4) + 
  #facet_grid(~facet, scales = "free", space = "free")+
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 10),
        axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 11, angle = 45),
        axis.text.y = ggtext::element_markdown(hjust = 1, vjust = 0.5, size = 11),
        strip.text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "bottom",
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("")+ylab("")





geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Group", y = "Mean Value ± SE") + 
  ggh4x::facet_grid2(cols = vars(module), rows = vars(rows), scales = "free_x", space = "free_x")+
  theme_classic()+ 
  #ggh4x::facet_grid2(cols = vars(module.to))+
  # scale_alpha_manual(values = c(0.5, 0.75, 1)) +
  scale_fill_gradient2(low = "grey100",  mid = "grey50", high = "red", midpoint = 0.4) +
  geom_text(aes(label = letters, y = mean + sd + 0.05), 
            angle = 0, hjust = 0.5, vjust = 0.5, size = 4, color = "black") +
  theme(panel.spacing = unit(2, "mm"),        
        axis.text.y.left  = element_text(size = 10),
        axis.text.x = ggtext::element_markdown(hjust = 1, vjust = 1, size = 11, angle = 45),
        axis.text.y = ggtext::element_markdown(hjust = 1, vjust = 0.5, size = 11),
        strip.text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  xlab("")+ylab("")
barchart.average.cor.ME.contigs.water

