#Script 2 - Summary Stats
#set working directory 
setwd("/Volumes/USB20FD/PhD Projects. Meetings. Protocols/Developmental Calf.pooling/Individual Calf/Code for Pub/16S/")

#load the packages and the data
source("Script1_LoadData.R", echo=TRUE)

##lets look at the number of reads per sample and the distribution
sample_sum_df <- data.frame(sum = sample_sums(data))  
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 7500) +
  ggtitle("Distribution of total ASVs within each samples") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) # looks alright, skews left

calf <- subset_samples(data, sample_type=="Calf" & SampleID!="Indiv16S_30") #remove sample 30
calf # 46 individual calves only
calf <- prune_taxa(taxa_sums(calf) > 0, calf)

min(sample_sums(calf)) #858904
max(sample_sums(calf)) #2106460
median(sample_sums(calf)) #1254559
mean(sample_sums(calf)) #1362710

#Take sample sums and test if there are any differences in the sequencing depths
sequencing_reads_poolsvsindivs<- ggplot(sequencing_reads, aes(x = sample_type, y= forward_count, fill= sample_type, colour = sample_type)) +
  theme_bw() +
  labs(title = "SEQUENCED READS PER SAMPLE TYPE", y= "ASVs per Sample", color= "Sample Type", fill = "Sample Type") +
  # facet_wrap(~, scales = "free", labeller = as_labeller(c("DNA_repool" = "DNA Repools",
  #                                                               "Calf" = "Calves"))) +
  geom_boxplot(alpha = 0.1, position = position_dodge2()) +
  geom_point(size = 3, shape = 18, position = position_dodge(width= 0.75)) +
  scale_fill_manual(values = upsetcolor) +
  scale_color_manual(values = upsetcolor) +
  expand_limits(y = 0) +
  scale_y_continuous(expand= c(0.0012,0,0.1,0),labels = scales::label_comma()) +
  scale_x_discrete(labels=c("Calf"="Calves", "DNA_pool"="DNA Pool","Fecal_pool"="Fecal Pool"))+
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust= 0.95),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.5),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc(method = "wilcox_test", p.adjust.method = "BH", label = "p = {p.adj}",
           hide.ns = TRUE,
           step.increase = 0.08)
sequencing_reads_poolsvsindivs

# #Do this next - keep in mind the y=sum needs to change; you would end up with a sum by sample type  
# sequencing_depth_poolsvsindivs<- ggplot(#put new md object from script 1 sample_sums_md$ASV_counts, aes(x = sample_type, y= sum, fill= sample_type, colour = sample_type)) +
#   theme_bw() +
#   labs(title = "CLASSIFIED READ DEPTH PER SAMPLE TYPE", y= "ASVs per Sample", color= "Sample Type", fill = "Sample Type") +
#   geom_boxplot(alpha = 0.1, position = position_dodge2()) +
#   geom_point(size = 3, shape = 18, position = position_dodge(width= 0.75)) +
#   scale_fill_manual(values = upsetcolor) +
#   scale_color_manual(values = upsetcolor) +
#   expand_limits(y = 0) +
#   scale_y_continuous(expand= c(0.0012,0,0.1,0)) +
#   scale_x_discrete(labels=c("Calf"="Calves", "DNA_pool"="DNA Pool","Fecal_pool"="Fecal Pool"))+
#   theme(legend.position = "none",
#         plot.title = element_text(size = 32),
#         panel.border = element_rect(colour = "black", linewidth= 1),
#         strip.background = element_rect(fill = "black"),
#         strip.text = element_text(colour = "white", size = 16, face = "bold"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust= 0.95),
#         axis.ticks.x = element_line(colour = "black", linewidth = 0.5),
#         axis.title.y = element_text(size = 32),
#         axis.text.y = element_text(colour = "black", size = 14),
#         axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
#   geom_pwc(method = "wilcox_test", p.adjust.method = "BH", label = "p = {p.adj}",
#            hide.ns = TRUE,
#            step.increase = 0.08)
# sequencing_depth_poolsvsindivs

########Agglomerate at different levels of classification#### Summary Stats 
#make a dataframe
taxacalf.df<-as.data.frame(tax_table(calf))

#filter unclassified phylum
phylum.ps <- tax_glom(calf, taxrank = "Phylum", NArm = F)
phylum.ps

unclassified_phylum.df<-taxacalf.df %>% filter(grepl('unclassified',Phylum))
unclassified_phylum.df
#Pull out unique
unclassified_phylum <- row.names(unclassified_phylum.df)
unclassified_phylum

# Example of a function that we can use (and re-use) to remove unwanted taxa
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

# We can run the "pop_taxa" function here using the phylum phyloseq object and the list of features we want to remove
trimmed_phylum.ps = pop_taxa(phylum.ps, unclassified_phylum)

sum(sample_sums(trimmed_phylum.ps))  /  sum(sample_sums(calf))*100

#calculate the % unclassified Phyla
total_reads <- sum(sample_sums(phylum.ps))

unclassified_reads <- sum(
  taxa_sums(prune_taxa(unclassified_phylum, phylum.ps))
)

percent_unclassified <- unclassified_reads / total_reads * 100
percent_unclassified #9.09 e-5

###filter unclassified classes
class.ps <- tax_glom(calf, taxrank = "Class", NArm = F)
class.ps

unclassified_class.df<-taxacalf.df %>% filter(grepl('unclassified',Class))
unclassified_class.df
#Pull out unique
unclassified_class <- row.names(unclassified_class.df)
unclassified_class
#pull out the unclassifieds
trimmed_class.ps = pop_taxa(class.ps, unclassified_class)

sum(sample_sums(trimmed_class.ps))  /  sum(sample_sums(calf)) * 100

###filter unclassified order
order.ps <- tax_glom(calf, taxrank = "Order", NArm = F)
order.ps
unclassified_order.df<-taxacalf.df %>% filter(grepl('unclassified',Order))
unclassified_order.df
#Pull out unique
unclassified_order <- row.names(unclassified_order.df)
unclassified_order
#pull out the unclassifieds
trimmed_order.ps = pop_taxa(order.ps, unclassified_order)

sum(sample_sums(trimmed_order.ps))  /  sum(sample_sums(calf)) * 100

###filter out unclassified family
family.ps <- tax_glom(calf, taxrank =  "Family", NArm=F)
family.ps
unclassified_family.df<-taxacalf.df %>% filter(grepl('unclassified',Family))
unclassified_family.df
#Pull out unique
unclassified_family <- row.names(unclassified_family.df)
unclassified_family
#pull out the unclassifieds
trimmed_family.ps = pop_taxa(family.ps, unclassified_family)

sum(sample_sums(trimmed_family.ps))  /  sum(sample_sums(calf)) * 100

###Filter out unclassified genus
genus.ps <- tax_glom(calf, taxrank = "Genus", NArm=F)
genus.ps
unclassified_genus.df<-taxacalf.df %>% filter(grepl('unclassified',Genus))
unclassified_genus.df
#Pull out unique
unclassified_genus <- row.names(unclassified_genus.df)
unclassified_genus
#pull out the unclassifieds
trimmed_genus.ps = pop_taxa(genus.ps, unclassified_genus)

sum(sample_sums(trimmed_genus.ps))  /  sum(sample_sums(calf)) * 100

###Filter unclassified species 
species.ps <- tax_glom(calf, taxrank = "Species", NArm=F)
species.ps
unclassified_species.df<-taxacalf.df %>% filter(grepl('unclassified',Species))
unclassified_species.df
#Pull out unique
unclassified_species <- row.names(unclassified_species.df)
unclassified_species
#pull out the unclassifieds
trimmed_species.ps = pop_taxa(species.ps, unclassified_species)

sum(sample_sums(trimmed_species.ps))  /  sum(sample_sums(calf)) * 100

# 16S Rarefaction Curve -------------------------------------------------------
#subset ASV table
otu_matrix <- otu_table(calf)
otu_matrix <- as.matrix(otu_matrix)


#Calculate rarefaction curve
set.seed(42)

calculate_rarefaction_curves <- function(otu_matrix, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(otu_matrix, measures, depth) {
    if(max(sample_sums(otu_matrix)) < depth) return()
    otu_matrix <- prune_samples(sample_sums(otu_matrix) >= depth, otu_matrix)
    
    rarified_otu_matrix <- rarefy_even_depth(otu_matrix, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_otu_matrix, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, otu_matrix = otu_matrix, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

#calculate observed and Shannon at different depths (caps at 1 mil)
rarefaction_curve_data <- calculate_rarefaction_curves(otu_matrix, c('Observed'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(calf)), by.x = 'Sample', by.y = 'SampleID')

#Okay back to your regularly scheduled programming
View(rarefaction_curve_data_summary_verbose)

# Create the paired plot showing the rarefaction curve data 
cohort_palette <- c("Pre3"="#559e83","Pre1"="#5a5255","Pre2"="#1b85b8","Post1"="#ae5a41","Post2"="#c3cb71")
rarefaction16s <-ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = weaning,
    group = Sample
  )
) +
  geom_line(size = 1) +
  scale_x_continuous(
    breaks = c(
      seq(0, 300000, by = 50000),
      seq(400000, 1000000, by = 100000)
    ),
    labels = comma
    #breaks = seq(0, max(amr_rarefaction_curve_data_summary_verbose$Depth), by = 100000),
    #labels = comma
  ) +
  facet_wrap(
    facets = ~ Measure,
    scales = 'free_y'
  ) +
  scale_colour_manual(values = cohort_palette) +  # Apply the custom color palette
  theme_bw() +
  labs(y = "Observed Richness", x = "Total ASVs Classified per Sample", title = "Rarefaction Curves by Cohort")+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 24, vjust = 1.5),
        axis.text.x = element_text(size = 24, colour = "black", angle=45, hjust=1),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),  # removes any remaining gray fill
        strip.background = element_blank(),   # removes gray bar
        strip.text = element_blank()
  )
rarefaction16s
ggsave("output_data/rarefactioncurve16s.png",rarefaction16s, width=)

# #Try a rarefaction curve looking at read depths - overlay the read depth onto the rarefaction curve 
# # Add a vertical line at each sample's original read depth (optional)
# read_depths <- data.frame(Sample = sample_names(otu_matrix),
#                           MaxDepth = sample_sums(otu_matrix))
# 
# ggplot(rarefaction_curve_data, aes(x = Depth, y = Alpha_diversity, color = Sample)) +
#   geom_line() +
#   geom_vline(data = read_depths, aes(xintercept = MaxDepth, color = Sample), linetype = "dashed")
# 
# #Try this chatgpt plot on just read depth 
# read_depths <- data.frame(Sample = sample_names(otu_matrix),
#                           ReadDepth = sample_sums(otu_matrix))
# 
# ggplot(read_depths, aes(x = Sample, y = ReadDepth)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_text(angle = 90))
# 
