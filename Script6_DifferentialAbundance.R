#Script for the Differential Abundance and Tornado plot
#ANCOMBC2 with class rows, group values
set.seed(123)

## Setup FECAL samples ####
# Specify the order of variables in a factor
# your first variable here will be used as the "reference" in the model
sample_data(calf)$weaning <- factor(sample_data(calf)$weaning ,levels = c("Post2", "Post1","Pre3","Pre2","Pre1"))

# Identify the taxa to keep by filtering out unwanted families
filtered_taxa <- taxa_names(calf)[!phyloseq::tax_table(calf)[, "Family"] %in% c("unclassified Bacteria", "unclassified Unassigned" ,"uncultured")]

# Prune the phyloseq object to keep only the selected taxa
calf <- prune_taxa(filtered_taxa, calf)

##### ANCOMBC - calf samples, comparing by weaning animal_num as fixed effect #####
###

# ANCOMBC-specific issue (not required for microbiome data):
# ANCOMBC looks for a "tax_level" to aggregate, and it expects it to match one of the typical taxonomic
# levels. This is not Pre3ortant if you aggregate the phyloseq object first.
# Change label for most specific annotation level, in this case micro gene group, to "Family"
#

# ANCOMBC model code
# Adjust the "fix_formula" to match your study design".
# Change the "group" variable to whatever you want to run pairwise comparisons on.
# We can play around with the "rand_formula" later for repeated measures.
ancom_output_calf = ancombc2(data = calf, assay_name = "counts", tax_level = "Family", 
                             fix_formula = "weaning", rand_formula = NULL,
                             p_adj_method = "holm", prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                             group = "weaning", struc_zero = TRUE, neg_lb = TRUE,
                             alpha = 0.05, n_cl = 2, verbose = TRUE,
                             global = TRUE, pairwise = TRUE,
                             dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
pairwise_ancom_output_calf <- ancom_output_calf$res_dunn

# Number of taxa in logFC comparisons
length(pairwise_ancom_output_calf$taxon) # 43


# Count the number of NAs per row in bias corrected taxa
# this shows in how many samples the taxa could not be adjusted for compositional effects
na_counts_per_row <- rowSums(is.na(ancom_output_calf$bias_correct_log_table))


setdiff(names(na_counts_per_row), pairwise_ancom_output_calf$taxon)

# Check the structure of the pairwise results
str(pairwise_ancom_output_calf)

# Filter the results to focus on comparisons between weanings
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_pairwise_ancom_output_calf = pairwise_ancom_output_calf %>%
  dplyr::select(taxon, contains("weaning")) 


df_fig_weaning_calf <- df_pairwise_ancom_output_calf %>%
  # Step 1: Apply the condition for the lfc_weaning columns
  dplyr::mutate(across(starts_with("lfc_weaning"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x,.x)),
                # Step 2: Round the lfc_weaning columns and create _rounded columns
                across(starts_with("lfc_weaning"), ~round(.x, 3), .names = "{.col}_rounded"),
                # Step 3: Assign colors based on diff_weaning and passed_ss columns
                across(starts_with("diff_weaning"), 
                       ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "darkred", "darkgrey"), 
                       .names = "{.col}_color")) %>%
  # Step 4: Pivot the _rounded columns into long format
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  # Step 5: Pivot the _color columns into long format
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  # Step 6: Filter to ensure the group and color_group match
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  # Step 7: Select the relevant columns and arrange by taxon
  select(-color_group) %>%
  arrange(taxon)

# Check which levels are available for comparison
as.factor(df_fig_weaning_calf$group)

df_fig_weaning_calf_filtered <- df_fig_weaning_calf %>%
  dplyr::filter(str_ends(group, "_rounded")) %>%
  dplyr::mutate(group = dplyr::case_when(
    group == "weaningPost1_rounded" ~ "Post1 vs Post2",
    group == "weaningPre3_rounded" ~ "Pre3 vs Post2",
    group == "weaningPre2_rounded" ~ "Pre2 vs Post2",
    group == "weaningPre1_rounded" ~ "Pre1 vs Post2",
    TRUE ~ group)) %>%
  dplyr::filter(group %in% c("Post1 vs Post2", 
                             "Pre3 vs Post2", 
                             "Pre2 vs Post2",
                             "Pre1 vs Post2")) %>%
  dplyr::mutate(group = factor(group, levels = c("Post1 vs Post2", 
                                                 "Pre3 vs Post2", 
                                                 "Pre2 vs Post2",
                                                 "Pre1 vs Post2"))) %>%
  droplevels()



# Confirm groups
levels(as.factor(df_fig_weaning_calf_filtered$group))

#### calf Dotplot - but with bias corrected abundances #######
# Assuming 'ancom_output_calf$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_calf <- ancom_output_calf$bias_correct_log_table

# Step 1: Replace NA values with the mean of the row
# Calculate row means, ignoring NAs
row_means <- rowMeans(log_table_calf, na.rm = TRUE)

# Replace NAs with corresponding row mean
log_table_calf <- as.data.frame(log_table_calf)  # Ensure it's a data frame for row-wise operations
for (i in seq_len(nrow(log_table_calf))) {
  log_table_calf[i, is.na(log_table_calf[i, ])] <- row_means[i]
}

sum(is.na(log_table_calf))

# Step 2: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_calf <- exp(log_table_calf)

# Step 3: Scale the pseudo-counts to a desired total count per sample
# Calculate the total sum of counts for each sample before scaling
total_counts_calf <- apply(pseudo_counts_calf, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_calf <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_calf <- t(t(pseudo_counts_calf) * (desired_total_count_calf / total_counts_calf))

# Convert scaled_counts_calf to a data frame
scaled_counts_calf.df <- as.data.frame(scaled_counts_calf)

# Add a column for taxa names
scaled_counts_calf.df$taxon <- rownames(scaled_counts_calf.df)

# Convert from wide to long format using pivot_longer
scaled_counts_calf.long <- scaled_counts_calf.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
calf.micro_data.melt <- psmelt(calf)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_calf.micro_data.melt <- calf.micro_data.melt %>%
  left_join(scaled_counts_calf.long, by = c("Family" = "taxon", "Sample" = "Sample"))

# Filter to include only rows with actual values (non-NA) in BiasAdj_counts
BiasAdj_calf.micro_data.filtered <- BiasAdj_calf.micro_data.melt %>%
  filter(!is.na(BiasAdj_counts))

# Count rows before and after filtering
before <- length(unique(BiasAdj_calf.micro_data.melt$Family))
after <- length(unique(BiasAdj_calf.micro_data.filtered$Family))

# Print the difference
cat("Family taxa removed:", before - after, "\n")
cat("Number of Family taxa modeled:", after, "\n")


Family_avg_counts_calf_by_weaning <- BiasAdj_calf.micro_data.filtered %>%
  dplyr::group_by(Family, weaning) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(weaning_label = dplyr::case_when(
    weaning == 'Post1' ~ "Post1 vs Post2",
    weaning == 'Pre3' ~ "Pre3 vs Post2",
    weaning == 'Pre2' ~ "Pre2 vs Post2",
    weaning == 'Pre1' ~ "Pre1 vs Post2",
    TRUE ~ as.character(weaning)  # Keep other weaning values unchanged
  ))


df_fig_weaning_calf_filtered_merged_data <- df_fig_weaning_calf_filtered %>%
  left_join(Family_avg_counts_calf_by_weaning, by = c("group" = "weaning_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_calf <- as.data.frame(phyloseq::tax_table(calf))

# Join the taxonomic data with your data frame
df_with_taxonomy_calf <- df_fig_weaning_calf_filtered_merged_data %>%
  left_join(tax_data_calf, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_calf <- df_with_taxonomy_calf %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_calf$group = factor(df_with_taxonomy_calf$group, levels = c("Post1 vs Post2", 
                                                                             "Pre3 vs Post2", 
                                                                             "Pre2 vs Post2",
                                                                             "Pre1 vs Post2"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_calf <- df_with_taxonomy_calf %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_calf <- ggplot(df_with_taxonomy_calf, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_calf %>% dplyr::filter(color == "darkred"),
                  aes(label = taxon),
                  color = "darkred",
                  size = 3,                  # Adjust label size as needed
                  position = position_dodge(width = 0.5),
                  box.padding = 0.4,         # Moderate padding to allow some vertical movement
                  point.padding = 0.4,       # Padding around points for some horizontal spacing
                  max.overlaps = Inf,        # Allow ggrepel to try placing all labels
                  segment.color = "grey70") +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = "Class", title = "logFC of Families") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size=15, face="bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black", face="bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 4, scales = "fixed")

main_species_plot_calf

# Create the taxonomy plot data
taxonomy_plot_data_calf <- df_with_taxonomy_calf %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_calf <- taxonomy_plot_data_calf %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_calf <- ggplot(taxonomy_plot_data_calf) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_calf$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_calf  <- plot_grid(
  taxonomy_plot_calf  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_calf  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.8,1.8),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_calf
# Save the final plot
ggsave("output_data//DA_dunnet_micro_family_calf_by_weaning.png", plot = combined_plot_calf, width = 20, height = 11, dpi = 300)
