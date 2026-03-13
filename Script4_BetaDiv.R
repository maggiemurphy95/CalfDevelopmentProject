#Script 4 - Beta-Diversity 
#############################################################################################
##############################         BETA-DIVERSITY         ##############################
#############################################################################################
#############################################################################################

calf.css <- phyloseq_transform_css(calf, log = F) # cumulative sum scaling normalization
calf.css

# check our read counts after normalization
plot(sort(sample_sums(calf.css), T), type = "h", ylab = "reads", xlab= "samples") #Skews left 

## calculate UniFrac dists
all_gunifrac.dist <- gunifrac(calf.css)

# ordinate
all_gunifrac.ord <- ordinate(calf.css, method = "NMDS", distance = all_gunifrac.dist)

# plot generalized
gen_ord <-plot_ordination(calf.css, all_gunifrac.ord, type = "sample", color = "weaning") + theme_bw() +
  labs(#title = "Microbiome",
    x= "NMDS1",
    y= "NMDS2") +
  geom_point(size = 3) +
  stat_ellipse(geom= "polygon", lty = 2, alpha = 0.18, aes(fill= weaning), size = 1, level = 0.99) +
  scale_colour_manual(values = cohort_palette) +
  scale_fill_manual(values = cohort_palette) +
  coord_fixed(ratio=2.5)+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 32),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(size = 0.9, colour = "black"))
gen_ord
ggsave("output_data/gen_ordination.png",gen_ord)

###############stats on overall general Ordination#######
beta_all.df <- as(sample_data(calf.css), "data.frame")
general.adonis.all <- adonis2(all_gunifrac.dist ~ weaning, beta_all.df)
general.adonis.all # singificant
general.pairwise_adonis.all <- pairwise.adonis2(all_gunifrac.dist ~ weaning, beta_all.df, perm = 9999, p.adjust.m = "BH")
general.pairwise_adonis.all # All significant but PW1 and PW2
general.disper <- betadisper(all_gunifrac.dist, beta_all.df$weaning)
general.permdisp <- permutest(general.disper, permutations = 9999, pairwise = T)
general.permdisp # NS;
