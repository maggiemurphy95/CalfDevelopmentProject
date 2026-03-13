#Script 5 - Dendrograms and Relative Abundance Plots 
######################dendrogram###
### Wards clustering on gunifrac
data.indivhclust <- hclust(all_gunifrac.dist, method = "ward.D2")
plot(data.indivhclust)
data.dendro.indiv <- as.dendrogram(data.indivhclust)
data.dendro.indivcalf <- dendro_data(data.dendro.indiv, type = "rectangle")
data.dendro.indivcalf

## add in metadatacolumn
write.csv(data.dendro.indivcalf[["labels"]][["label"]], "output_data/dendro_sample_order.csv")

weaning_col <- data.frame(weaning =c("Pre1",	"Pre1",	"Pre1",	"Pre1",	"Pre1",	"Pre1",	"Pre1",	"Post2",	"Pre3",	"Post1",	"Post1",	"Post2",	"Pre3",	"Post1",	"Post2",	"Post2",	"Post1",	"Post1",	"Post1",	"Post1",	"Post1",	"Post1",	"Post1",	"Post2",	"Post2",	"Post1",	"Post2",	"Post2",	"Post1",	"Pre3",	"Pre3",	"Pre3",	"Pre3",	"Pre3",	"Pre3",	"Pre3",	"Pre2",	"Pre2",	"Pre2",	"Pre2",	"Pre2",	"Pre2",	"Pre2",	"Pre2",	"Pre2",	"Pre2"))
age_col <- data.frame(age =c("3",	"3",	"2",	"2",	"3",	"2",	"2",	"94",	"90",	"92",	"88",	"97",	"90",	"90",	"95",	"95",	"92",	"93",	"90",	"90",	"90",	"89",	"93",	"99",	"99",	"93",	"97",	"97",	"87",	"90",	"93",	"90",	"93",	"90",	"90",	"93",	"37",	"36",	"39",	"37",	"36",	"39",	"36",	"39",	"39",	"39"))
data.dendro.indivcalf$labels <- cbind(data.dendro.indivcalf$labels, weaning_col, age_col)
data.dendro.indivcalf

#Plot the dendrogram
ward_dendro<-ggplot(data.dendro.indivcalf$segments) + theme_bw() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.5, lineend = "round", linejoin = "round") +
  geom_point(data = data.dendro.indivcalf$labels, aes(x,y, color = weaning), 
             size =9, position = position_nudge(y=-.08), shape = 15) +
  geom_text(data = data.dendro.indivcalf$labels,
            aes(x, y, label = age),  # replace label_variable with your actual label variable
            size = 4, position = position_nudge(y = -0.08)) +
  scale_colour_manual(values = cohort_palette) +
  theme(#legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(size = 0.75),
    axis.title.y = element_text(size = 36),
    axis.text.y = element_text(size = 20, colour = "black"))
ward_dendro
ggsave("output_data/ward_endro.png",ward_dendro)

###RA bar plot for under dendrogram
#Create an object for later use in all RA plots 
rel_abun_calf <- transform_sample_counts(calf.css,function(x) {x/sum(x)}*100)
# check our read counts after normalization, double check- should all be 100%
plot(sort(sample_sums(rel_abun_calf), TRUE), type = "p", ylab = "reads", xlab= "samples")

##Dendrogram specific RA plot 
dendrosampleorder<- c("Indiv16S_09",	"Indiv16S_08",	"Indiv16S_03",	"Indiv16S_05",	"Indiv16S_10",	"Indiv16S_04",	"Indiv16S_06",	"Indiv16S_46",	"Indiv16S_23",	"Indiv16S_34",	"Indiv16S_44",	"Indiv16S_42",	"Indiv16S_26",	"Indiv16S_31",	"Indiv16S_43",	"Indiv16S_48",	"Indiv16S_35",	"Indiv16S_33",	"Indiv16S_37",	"Indiv16S_38",	"Indiv16S_40",	"Indiv16S_36",	"Indiv16S_32",	"Indiv16S_49",	"Indiv16S_45",	"Indiv16S_50",	"Indiv16S_41",	"Indiv16S_47",	"Indiv16S_39",	"Indiv16S_25",	"Indiv16S_28",	"Indiv16S_22",	"Indiv16S_27",	"Indiv16S_21",	"Indiv16S_24",	"Indiv16S_29",	"Indiv16S_15",	"Indiv16S_17",	"Indiv16S_12",	"Indiv16S_16",	"Indiv16S_19",	"Indiv16S_13",	"Indiv16S_18",	"Indiv16S_20",	"Indiv16S_11",	"Indiv16S_14")
Ordercolors <- distinctColorPalette(154)

#Agglomerate and melt at Order 
RACalf_Order <- tax_glom(rel_abun_calf, taxrank = "Order", NArm = F)
RACalf_Order_melt <- psmelt(RACalf_Order)

dendro_ra<-ggplot(RACalf_Order_melt, aes(x= Sample, y= Abundance, fill = Order)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = Ordercolors) +
  scale_x_discrete(limits = dendrosampleorder) +
  theme(legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.75),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 0.75, lineend = "square", colour = "black"),
    axis.title.y = element_text(size = 36),
    axis.text.y = element_text(size = 20, colour = "black"))
dendro_ra
ggsave("output_data/raplot_dendro.png",dendro_ra)
#Use BioRender to combine the two plots

#Average calculations on Order RA
#Use the rel_abun_calf objects
# calculate percentages across all order counts for calves
calf_order_abund_perc <- RACalf_Order_melt %>%
  dplyr::group_by(weaning, Order) %>%
  dplyr::summarize(mean_order_calf_perc = mean(Abundance), se_order_calf_perc = (sd(Abundance)/sqrt(length(Abundance)))) %>%
  arrange(-mean_order_calf_perc)
write.csv(calf_order_abund_perc, "output_data/calf_order_abund_perc.csv")

######################Other Relative Abundances and Plots #######
#Using the rel_abun_calf object from above

#To make these Family plots, start by making a Phyla object
RACalf_Phyla <- tax_glom(rel_abun_calf, taxrank = "Phylum", NArm = F)
RACalf_Phyla_melt <- psmelt(RACalf_Phyla)

#Make a Family object
RACalf_Family <- tax_glom(rel_abun_calf, taxrank = "Family", NArm = F)
RACalf_Family_melt <- psmelt(RACalf_Family)

########################################Firmicutes Families Graph################################################
FirmicutesRA <- subset_taxa(RACalf_Phyla, Phylum=="Firmicutes") %>%
  psmelt()

#Firmicutes Families graph from Lee phyloseq script 
firmicutes_family_calf <- subset_taxa(RACalf_Family, Phylum=="Firmicutes") %>%
  psmelt()

firmicutes_phylum_calf <- subset_taxa(RACalf_Phyla, Phylum=="Firmicutes") %>%
  psmelt()

firmicutes_family_palette <- distinctColorPalette(70)

firmicutes_family_calf <- ggplot(firmicutes_family_calf, aes(x= weaning, y= Abundance)) +
  theme_bw() + coord_flip() +
  labs(y= "Relative Abundance (%)", title = "Firmicutes") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black", size = 0.75) +
  geom_errorbar(firmicutes_phylum_calf,
                mapping = aes(x= weaning, y= Abundance), 
                stat = "summary", size = 0.75, width = 0.4) +
  scale_fill_manual(values = firmicutes_family_palette) +
  scale_y_continuous(limits = c(0,75), expand = c(0.001,0,-0.1,0)) +
  scale_x_discrete(limits = c("Post2","Post1","Pre3","Pre2","Pre1")) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 36, vjust = -1),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 32),
        axis.ticks.x = element_line(size = 0.9, colour = "black"),
        axis.ticks.y = element_line(size = 0.9, colour = "black"))=
firmicutes_family_calf
ggsave("output_data/firmicutes_family.png",firmicutes_family_calf)

#Percentages 
# calculate percentages across firmicutes family counts for calves
calf_abund_perc_firm <- firmicutes_family_calf %>%
  dplyr::group_by(weaning, Family) %>%
  dplyr::summarize(mean_order_calf_perc = mean(Abundance), se_order_calf_perc = (sd(Abundance)/sqrt(length(Abundance)))) %>%
  arrange(-mean_order_calf_perc)
write.csv(calf_abund_perc_firm, "output_data/calf_FirmicutesFAM_abund_perc.csv")


########################################Proteobacteria Families Graph################################################
proteobacteria_family_calf <- subset_taxa(RACalf_Family, Phylum=="Proteobacteria") %>%
  psmelt()
proteobacteria_phylum_calf <- subset_taxa(RACalf_Phyla, Phylum=="Proteobacteria") %>%
  psmelt()
proteobacteria_family_palette <- distinctColorPalette(65)

proteo_family_calf <-ggplot(proteobacteria_family_calf, aes(x= weaning, y= Abundance)) +
  theme_bw() + coord_flip() +
  labs(y= "Relative Abundance (%)", title = "Proteobacteria") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black", size = 0.75) +
  geom_errorbar(proteobacteria_phylum_calf,
                mapping = aes(x= weaning, y= Abundance), 
                stat = "summary", size = 0.75, width = 0.4) +
  scale_fill_manual(values = proteobacteria_family_palette) +
  scale_x_discrete(limits = c("Post2", "Post1","Pre3", "Pre2","Pre1"), labels = c("Post2", "Post1","Pre3", "Pre2","Pre1")) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 36, vjust = -1),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 32),
        axis.ticks.x = element_line(size = 0.9, colour = "black"),
        axis.ticks.y = element_line(size = 0.9, colour = "black"))
proteo_family_calf
ggsave("output_data/proteo_family_calf.png",proteo_family_calf)

#Percentages 
# calculate percentages across proteobacteria order counts for calves
calf_abund_perc_proteo <- proteobacteria_family_calf %>%
  dplyr::group_by(weaning, Family) %>%
  dplyr::summarize(mean_order_calf_perc = mean(Abundance), se_order_calf_perc = (sd(Abundance)/sqrt(length(Abundance)))) %>%
  arrange(-mean_order_calf_perc)
write.csv(calf_abund_perc_proteo, "output_data/calf_ProteoFAM_abund_perc.csv")

##########################################Bacteroidetes Families Graph#################################################
bacteroidetes_family_calf <- subset_taxa(RACalf_Family, Phylum=="Bacteroidota") %>%
  psmelt()
bacteroidetes_phylum_calf <- subset_taxa(RACalf_Phyla, Phylum=="Bacteroidota") %>%
  psmelt()
bacteroidetes_family_palette <- distinctColorPalette(53)

bacteroid_family_calf<-ggplot(bacteroidetes_family_calf, aes(x= weaning, y= Abundance)) +
  theme_bw() + coord_flip() +
  labs(y= "Relative Abundance (%)", title = "Bacteroidota") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black", size = 0.75) +
  geom_errorbar(bacteroidetes_phylum_calf,
                mapping = aes(x= weaning, y= Abundance), 
                stat = "summary", size = 0.75, width = 0.4) +
  scale_fill_manual(values = bacteroidetes_family_palette) +
  scale_y_continuous(limits = c(0,75), expand = c(0.001,0,-0.1,0)) +
  scale_x_discrete(limits = c("Post2","Post1","Pre3","Pre2","Pre1"), labels = c("Post2","Post1","Pre3","Pre2","Pre1")) +
  theme(legend.position = "none",
    plot.margin = unit(c(1,1,1,1),"lines"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 40),
    axis.title.x = element_text(size = 36, vjust = -1),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 32),
    axis.ticks.x = element_line(size = 0.9, colour = "black"),
    axis.ticks.y = element_line(size = 0.9, colour = "black"))
bacteroid_family_calf
ggsave("output_data/bacteroid_family_calf.png",bacteroid_family_calf)

#Percentages 
# calculate percentages across all order counts for calves
calf_abund_perc_bacteroid <- bacteroidetes_family_calf %>%
  dplyr::group_by(weaning, Family) %>%
  dplyr::summarize(mean_order_calf_perc = mean(Abundance), se_order_calf_perc = (sd(Abundance)/sqrt(length(Abundance)))) %>%
  arrange(-mean_order_calf_perc)
write.csv(calf_abund_perc_bacteroid, "output_data/calf_BacteroidFAM_abund_perc.csv")