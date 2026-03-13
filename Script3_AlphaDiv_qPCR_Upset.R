#Script 3 - Alpha-Diversity; qPCR data and UpsetPlots
#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

## data preparation - calculation and adding metadata
alpha_div1 <- estimate_richness(calf, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2 <- estimate_pd(calf) # calculating Faith's PD

alpha_div <- cbind(alpha_div1, alpha_div2) # combining Faith's and other metrics
alpha_div # looks good, but have some duplicates so lets trim it a bit
alpha_div <- alpha_div[,c(1:5)]
alpha_div # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div.df <- as(sample_data(calf), "data.frame")

alpha_div_calf <- cbind(alpha_div, alpha_div.df)

#Richness
#box plot
richness16s_boxplot<-ggplot(alpha_div_calf, aes(x=weaning, y= Observed, fill= weaning, color=weaning)) + 
  theme_bw() +
  geom_boxplot(alpha= 0.5, size= 0.75)+
  geom_point()+
  labs(y= "Richness", x= "") +
  scale_fill_manual(values = cohort_palette)+
  scale_color_manual(values = cohort_palette)+
  scale_x_discrete(limits=c("Pre1","Pre2","Pre3","Post1","Post2")) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
richness16s_boxplot
ggsave("output_data/richness16s_boxplot.png",richness16s_boxplot)

#Pairwise Wilcoxon on alpha-div
pairwise.wilcox.test(alpha_div_calf$Observed, alpha_div_calf$weaning, p.adjust.method = "BH")

#Shannon/Diversity 
##Boxplot
diversity16s_boxplot<-ggplot(alpha_div_calf, aes(x=weaning, y= Shannon, fill= weaning, color=weaning)) + 
  theme_bw() +
  geom_boxplot(alpha= 0.5, size= 0.75)+
  geom_point()+
  labs(y= "Shannon's Diversity", x= "") +
  scale_fill_manual(values = cohort_palette)+
  scale_color_manual(values = cohort_palette)+
  scale_x_discrete(limits=c("Pre1","Pre2","Pre3","Post1","Post2")) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
diversity16s_boxplot
ggsave("output_data/diversity16s_boxplot.png",diversity16s_boxplot)

#Pairwise Wilcoxon on 16S Diversity
pairwise.wilcox.test(alpha_div_calf$Shannon, alpha_div_calf$weaning, p.adjust.method = "BH")

###qPCR DATA##########
#create the csv with the results will only work on single sheets 
qPCRresults <- read.csv("data/NAHMSqPCRmathredo.csv")

#shapirowilk test - Make sure your excel file has your means listed as numbers
shapiro.test(qPCRresults$newcopy) #4.82 e-06 so this is assumed to be not-normal

#compute analysis of variance for nonparametric data
kruskal.test(newcopy ~ weaning, data = qPCRresults)#0.3611 so no significance
#make it an object to use

qPCR_boxplot<-ggplot(qPCRresults, aes(x=weaning, y= newcopy, fill= weaning, color=weaning)) + 
  theme_bw() +
  geom_boxplot(alpha= 0.5, size= 0.75)+
  geom_point()+
  labs(title="Total Microbial Abundance", y= "Copies/mg feces", x= "") +
  scale_fill_manual(values = cohort_palette)+
  scale_color_manual(values = cohort_palette)+
  scale_x_discrete(limits=c("Pre1","Pre2","Pre3","Post1","Post2")) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
qPCR_boxplot
ggsave("output_data/qPCR_boxplot.png",qPCR_boxplot)

#####################UPSET PLOTS ON SHARED ASVs##############
#load MicrobiotaProcess now, but it conflicts with Phyloseq, so detach after. 
library("MicrobiotaProcess")
library(UpSetR)
#Do some filtering and set colors
#pre-filtering before making the upset plot - get rid of singletons and doubletons
data.preDA <- preDA(calf, min.samples = 2) #5939 grouped as others in output
upset_DA <- get_upset(data.preDA, factorNames="weaning")

#make a way to save the plot
png("output_data/upsetplot.png", width = 1200, height = 800)
#Then make an upset plot
upset(upset_DA, sets=c("Post2","Post1","Pre3","Pre2","Pre1"
),
sets.bar.color = c("goldenrod2","#ae5a41", "#559e83","#1b85b8","honeydew4"),
keep.order=TRUE, empty.intersections = "on",
order.by = "freq",
text.scale = c(3, 3, 2, 2, 3, 3))

dev.off()

#Now detach so not to conflict with Phyloseq later
detach("package:MicrobiotaProcess", unload = TRUE)
detach("package:UpSetR", unload = TRUE)#detach this package now so it doesn't impact phyloseq later
