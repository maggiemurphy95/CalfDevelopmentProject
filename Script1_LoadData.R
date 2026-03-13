#Script 1 - Load packages, source functions, load data 
#load packages
library(btools)
library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(metagenomeSeq)
library(metagMisc)
library(pairwiseAdonis)
library(phyloseq)
library(picante)
library(randomcoloR)
library(scales)
library(stringr)
library(vegan)
library(microbiome)
library(devtools)
library(ggdendro)
#library("MicrobiotaProcess") #Only load when doing UpsetPlots - conflicts with Phyloseq
library(UpSetR)
library(DAtest)
library(ANCOMBC)
library(DT)
library(ggrepel)
library(tidyr)
library(ggpubr)

#set working directory 
setwd("/Volumes/USB20FD/PhD Projects. Meetings. Protocols/Developmental Calf.pooling/Individual Calf/Code for Pub/16S/")

# source custom MEG lab scripts####
source("custommeg_scripts/g_unifrac.R")
source("custommeg_scripts/changeGGTaxaNames.R")

# import data####
qiimedata <- import_biom("data/table-with-taxonomy.biom", "data/tree.nwk", "data/dna-sequences.fasta")
#MM: Make sure all of your titles are lowercase
map_file <- import_qiime_sample_data("data/metadatafilenew.txt")
str(map_file)

# combining sample data with the rest
data <- merge_phyloseq(qiimedata,map_file)
#To see your variables and to make sure everything looks okay:
sample_variables(data)
## DATA EXPLORATION
data # we have 77 samples, but one failed because it only had 180 reads when normal is 135,000: Check DNAFecalPool13

# check the names of our ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data) # beauty, now they are named properly

# changing the GG style naming (k__Bacteria, etc.); Would not let me change- said it did not know str replace function
tax.data <- data.frame(tax_table(data)) # extract the taxonomy table as a data frame
tax.data.names <- changeGGtaxa(tax.data) # this gets rid of the GG format

# now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- kingdom
  } else if (tax.data.names[i,3] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- phylum
  } else if (tax.data.names[i,4] == ""){
    class <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- class
  } else if (tax.data.names[i,5] == ""){
    order <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- order
  } else if (tax.data.names[i,6] == ""){
    family <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- family
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Species[i] <- paste("unclassified ",tax.data.names$Genus[i], sep = "")
  }
}

head(tax.data.names) # great, no more NAs and no more k__
tax_table(data) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
tail(tax_table(data), 20) # sweet, lookin good!

# Assuming final object is "data"
sample_data(data)$ASV_counts <- sample_sums(data)

# in case you run into issues with extracting the metadata use this command
metadata.df <- as(sample_data(data), "data.frame")



