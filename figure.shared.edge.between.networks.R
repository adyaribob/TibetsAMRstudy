library(igraph)
library(dplyr)
library(igraph)

#### bac and euk ####
dat1 <- readRDS("dat.res.net.influent.rds")
dat1 <- dat1[!is.na(dat1$pval), ]
network.influent <- igraph::graph_from_data_frame(dat1[dat1$pval < 0.05 & dat1$Weight > 0, ], directed = FALSE, vertices = NULL)
network.influent <- simplify(network.influent)
a <- as_data_frame(network.influent, what = "edge")
a$id <- paste(a$from, a$to, sep = "___")
a <- merge(a, dat1)
a <- a[,c("Source","Target","Weight","pval","padj","from.group","to.group","relationship","id")]
network.influent <- igraph::graph_from_data_frame(a[a$relationship %in% c("euk_bac","bac_euk","bac_bac","euk_euk"), ], directed = FALSE)

dat1 <- readRDS("dat.res.net.water.rds")
dat1 <- dat1[!is.na(dat1$pval), ]
dat1$relationship %>% unique
network.water <- igraph::graph_from_data_frame(dat1[dat1$pval < 0.05 & dat1$Weight > 0, ], directed = FALSE, vertices = NULL)
network.water <- simplify(network.water)
a <- igraph::as_data_frame(network.water, what = "edge")

a$relationship %>% unique
a$id <- paste(a$from, a$to, sep = "___")
a <- merge(x = a, y = dat1, all.x = "TRUE", by = "id")
a <- a[,c("Source","Target","Weight","pval","padj","from.group","to.group","relationship","id")]
network.water <- igraph::graph_from_data_frame(a[a$relationship %in% c("euk_bac","bac_euk","bac_bac","euk_euk"), ], directed = FALSE)

dat1 <- readRDS("dat.res.net.sediment.rds")
dat1 <- dat1[!is.na(dat1$pval), ]
network.sediment <- igraph::graph_from_data_frame(dat1[dat1$pval < 0.05 & dat1$Weight > 0, ], directed = FALSE, vertices = NULL)
network.sediment <- simplify(network.sediment)
a <- igraph::as_data_frame(network.sediment, what = "edge")
a$id <- paste(a$from, a$to, sep = "___")
a <- merge(a, dat1)
a <- a[,c("Source","Target","Weight","pval","padj","from.group","to.group","relationship","id")]
network.sediment <- igraph::graph_from_data_frame(a[a$relationship%in% c("euk_bac", "bac_euk","bac_bac","euk_euk"), ], directed = FALSE)

#### water module and influent #####
a <- dat.modularity.water.bac.euk.waktrap
mod1<- intersection(network.influent, 
             subgraph(graph = network.water, 
                      vids = a[a$wtc.membership == 1, "wtc.names"])) %>% E %>% length

mod2<- intersection(network.influent, 
             subgraph(graph = network.water, 
                      vids = a[a$wtc.membership == 2, "wtc.names"])) %>% E %>% length

mod3<- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 3, "wtc.names"])) %>% E %>% length

mod4<- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 4, "wtc.names"])) %>% E %>% length

mod5<- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 5, "wtc.names"])) %>% E %>% length

mod6<- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 6, "wtc.names"])) %>% E %>% length

mod7<- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 7, "wtc.names"])) %>% E %>% length

mod8<- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 8, "wtc.names"])) %>% E %>% length
mod9 <- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 9, "wtc.names"])) %>% E %>% length

mod10 <- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 10, "wtc.names"])) %>% E %>% length

mod11 <- intersection(network.influent, 
                    subgraph(graph = network.water, 
                             vids = a[a$wtc.membership == 11, "wtc.names"])) %>% E %>% length


dat.shared.network.water.inf.baceuk <- data.frame(water.with.influent = paste("mod", 1:11, sep = ""), same.edge = c(mod1, mod2, mod3, mod4, mod5,
                                                                                                                    mod6, mod7, mod8, mod9, mod10, mod11))

dat.shared.network.water.inf.baceuk$tot.edge <- c(length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 1, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 2, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 3, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 4, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 5, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 6, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 7, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 8, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 9, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 10, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.water, 
                                                                    vids = a[a$wtc.membership == 11, "wtc.names"])))
                                                  )

#### sediment module and influent ####
dat.modularity.sediment.bac.euk.waktrap$wtc.membership %>% unique
a <- dat.modularity.sediment.bac.euk.waktrap
mod1<- intersection(network.influent, 
                    subgraph(graph = network.sediment, 
                             vids = a[a$wtc.membership == 1, "wtc.names"])) %>% E %>% length

mod2<- intersection(network.influent, 
                    subgraph(graph = network.sediment, 
                             vids = a[a$wtc.membership == 2, "wtc.names"])) %>% E %>% length

mod3<- intersection(network.influent, 
                    subgraph(graph = network.sediment, 
                             vids = a[a$wtc.membership == 3, "wtc.names"])) %>% E %>% length

mod4<- intersection(network.influent, 
                    subgraph(graph = network.sediment, 
                             vids = a[a$wtc.membership == 4, "wtc.names"])) %>% E %>% length

mod5<- intersection(network.influent, 
                    subgraph(graph = network.sediment, 
                             vids = a[a$wtc.membership == 5, "wtc.names"])) %>% E %>% length

dat.shared.network.sediment.inf.baceuk <- data.frame(sediment.with.influent = paste("mod", 1:5, sep = ""), 
           same.edge = c(mod1, mod2, mod3, mod4, mod5))


dat.shared.network.sediment.inf.baceuk$tot.edge <- c(length(E(subgraph(graph = network.sediment, 
                                                                    vids = a[a$wtc.membership == 1, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.sediment, 
                                                                    vids = a[a$wtc.membership == 2, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.sediment, 
                                                                    vids = a[a$wtc.membership == 3, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.sediment, 
                                                                    vids = a[a$wtc.membership == 4, "wtc.names"]))),
                                                  length(E(subgraph(graph = network.sediment, 
                                                                    vids = a[a$wtc.membership == 5, "wtc.names"])))
                                                  )
dat.shared.network.sediment.inf.baceuk$shared.per <- dat.shared.network.sediment.inf.baceuk$same.edge/dat.shared.network.sediment.inf.baceuk$tot.edge

#### barchart figure 3 ######

colnames(dat.shared.network.sediment.inf.baceuk)[1] <- "module"
colnames(dat.shared.network.water.inf.baceuk)[1] <- "module"
dat.shared.network.sediment.inf.baceuk$media <- "Sediment"
dat.shared.network.water.inf.baceuk$media <- "Water"

dat.shared.network.sediment.inf.baceuk$module[dat.shared.network.sediment.inf.baceuk$module == "mod1"] <- "Module 1"
dat.shared.network.sediment.inf.baceuk$module[dat.shared.network.sediment.inf.baceuk$module == "mod2"] <- "Module 2"
dat.shared.network.sediment.inf.baceuk$module[dat.shared.network.sediment.inf.baceuk$module == "mod3"] <- "Module 3"
dat.shared.network.sediment.inf.baceuk$module[dat.shared.network.sediment.inf.baceuk$module == "mod4"] <- "Module 4"
dat.shared.network.sediment.inf.baceuk$module[dat.shared.network.sediment.inf.baceuk$module == "mod5"] <- "Module 5"

dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod1"] <- "Module 1"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod2"] <- "Module 2"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod3"] <- "Module 3"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod4"] <- "Module 4"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod5"] <- "Module 5"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod6"] <- "Module 6"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod7"] <- "Module 7"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod8"] <- "Module 8"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod9"] <- "Module 9"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod10"] <- "Module 10"
dat.shared.network.water.inf.baceuk$module[dat.shared.network.water.inf.baceuk$module == "mod11"] <- "Module 11"


c <- rbind(dat.shared.network.sediment.inf.baceuk, dat.shared.network.water.inf.baceuk)
c$media <- factor(c$media, levels = c("Water","Sediment"))
c$module <- factor(c$module, levels = c("Module 1", "Module 2", "Module 3", "Module 4", "Module 5",
                                        "Module 6", "Module 7", "Module 8", "Module 9", "Module 10", "Module 11"))
dat.shared.network.sediment.water.with.inf.baceuk <- c
fig.shared.edge.water.sediment.with.influent.barchart <- ggplot(data = dat.shared.network.sediment.water.with.inf.baceuk,
                                                                aes(x = module, y = rel)) + 
  geom_bar(stat = "identity", aes(fill = rel)) +
  scale_fill_gradient2(low = "white", high = "#c88b5cff", name = "Similarity to sewage")+
  ggh4x::facet_grid2(cols = vars(media), scales = "free", space = "free") + 
  theme_classic()+labs(x = "", y = "Number of shared edges\nwith sewage network (%)")+
  theme(
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text.y.left  = element_text(size = 10, colour = "black"),
    strip.text.y = element_text(size = 9, colour = "black"), 
    axis.title.y = element_text(size = 10),
    legend.key = NULL,
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 10.5 ,face = "plain", hjust = 1, vjust = 1, angle = 45)
  )

