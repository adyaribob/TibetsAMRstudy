
dat.figS6.neutral.model <- plot.neutral.contig$data
dat.figS6.neutral.model.R2 <- plot.neutral.contig$layers[[5]]$data


plot.neutral.contig <- ggplot(dat.figS6.neutral.model, aes(x=log10(p), y=freq))+
  geom_point(aes(color = lefse),alpha=0.5, size=1.5)+
  geom_line(aes(x=log10(p), y=freq.pred),
            color="blue",linetype="solid",size=0.7)+
  geom_line(aes(x=log10(p), y=pred.lwr),
            color="blue",linetype="dashed",size=0.7)+
  geom_line(aes(x=log10(p), y=pred.upr),
            color="blue",linetype="dashed",size=0.7)+
  scale_color_manual(name = "LEfSe", values = c("SD" = "#747070ff",  "WA" = "#31b2e6ff",
                                                "IF" = "#c88b5cff" , "None" = "purple1"))+ 
  ggh4x::facet_grid2(cols = vars(media))+
  geom_text(data = dat.figS6.neutral.model.R2, 
            aes(x = p,
                y = freq,
                label = R),
            parse = TRUE, # Parse mathematical expressions
            #inherit.aes = FALSE,
            size = 4,
            color = "black")+
  xlab("Log10 (Mean Relative abundance)") + ylab("Frequency")+ theme_classic()+
  theme(
    strip.text.x = element_text(size = 12, colour = "black"), 
    strip.text.y = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    strip.background = element_rect(size = 0.5),
    strip.text.x.top  = element_text(size = 12),
    strip.text.y.right =   element_text(size = 12),
    legend.key = NULL,
    legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(size = 12 ,
                                           face = "plain", hjust = 1, vjust = 0.4))


figS6.neutral.contigs.barchart <- 
  ggplot(dat.figS6.neutral.model.barchart, aes(x = lefse, y = n, fill = upper_neutral_under)) +
  geom_bar(stat = "identity", position = position_stack()) +
  ylab("Number of contigs") + xlab("") + 
  scale_fill_manual(name = "lefse", values = c("above" = "#A52A2A",  "neutral" = "grey",
                                               "below" = "#29A6A6" ))+
  theme_classic()+   ggh4x::facet_grid2(pathogen ~ media, scales = "free_y", independent = "y")+
  scale_y_continuous(
    #limits = c(0, 100.05),          # Set y-axis range
    #breaks = seq(0, 100.05, 20),
    expand = c(0.01, 0)# Set labels at intervals of 0.2
  )+   
  geom_text(data = dat.figS6.neutral.model.barchart, 
            aes(label = n), position = position_stack(vjust = 0.5), size = 3) +  # Add text labels
  #  scale_fill_manual(values = c("sediment" = "#747070ff",  "water" = "#31b2e6ff",
  #                               "sewage" = "#c88b5cff", "none" = "grey80"), )+
  xlab("")+ 
  #stat_compare_means(label = "p.signif", size = 5, hide.ns = TRUE,  vjust = 0.4, 
  #                   comparisons = pairs )+
  theme(ggh4x.facet.nestline = element_line(linetype = 3),
        axis.text.y.left  = element_text(size = 12, colour = "black"),
        strip.text.x = element_text(size = 12, colour = "black"), 
        strip.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y = element_text(size = 12),
        strip.background = element_rect(size = 0.5),
        legend.key = NULL,
        legend.position = "bottom",
        axis.text.x = ggtext::element_markdown(size = 12 ,face = "plain", hjust = 1, vjust = 1, angle = 45))

