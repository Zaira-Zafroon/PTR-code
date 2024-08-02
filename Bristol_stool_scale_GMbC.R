####library####
library(dplyr)
library(purrr)
library(ggtext)
library(patchwork)
library(tidyverse)
library(rstatix)
######Cleaning and Importing#####
#Import the merged genome
merged_genome<-read.table("merged_genome.tsv",
                          sep = "\t", 
                          header = T,  comment.char = "",
                          check.names = F, na.strings = "NA")
#Import COPTR results
GMbC_2500<-read.csv("/Users/bioinfo/Documents/Masters/GMbC_2500_out.csv",
                    sep = ",", 
                    header = TRUE, 
                    comment.char = "",
                    check.names = FALSE, 
                    na.strings = "NA")
#head(GMbC_2500)

#colnames(GMbC_2500)

#Obtain the sepcies name from merged genome
subset_genome <- subset(merged_genome, select = c(user_genome,species))
#Change the column names for merging
names(GMbC_2500)[names(GMbC_2500) == "log2(PTR):genome_id/sample_id"] <- "genome_id"
names(subset_genome)[names(subset_genome) == "user_genome"] <- "genome_id"
#Clean the sample names
names(GMbC_2500) <- gsub("_R1_clean.fastq", "", names(GMbC_2500))
GMbC_2500 <- GMbC_2500 %>%
        left_join(subset_genome, by = "genome_id")
#unique(subset_genome$species)#4335
#Remove the duplicates
GMbC_2500 <- GMbC_2500 %>%
        distinct()
#Remove the Genome ID column
GMbC_2500 <- GMbC_2500 %>%
        select(-genome_id)
#Pivot Longer
GMbC_2500_long <- GMbC_2500 %>%
        pivot_longer(cols = -species, names_to = "SampleID", values_to = "PTR_value")
length(unique(GMbC_2500_long$species))#1942 species
#Remove NA 
GMbC_2500_long_no_na <- na.omit(GMbC_2500_long, cols = "PTR_value")
length(unique(GMbC_2500_long_no_na$species))

#Import the metadata file
GMbC_metadata<-read.table("metadata_gmbc_bn10_complete.tsv", sep = "\t", header = T, 
                          comment.char = "",check.names = F, na.strings = "NA")
head(GMbC_metadata)
#Create Sample Column
GMbC_metadata<- tibble::rownames_to_column(GMbC_metadata, "SampleID")
#Create a merged PTR and metdata file
GMbC_2500_long_no_na<-merge(GMbC_2500_long_no_na,GMbC_metadata,by="SampleID")
#Get the top 30 taxa
ranked_taxa<- read.table("ranked_taxa.tsv", sep = "\t", 
                         header = T,  comment.char = "",
                         check.names = F, na.strings = "NA")
top_30_taxa <- ranked_taxa %>% 
        slice(1:30)
# Renaming the column in the data frame
top_30_taxa  <- top_30_taxa  %>%
        rename(species = taxon)
####Bristol stool scale and relationship####
# Drop rows with missing Bristol Stool Scale values
GMbC_2500_long_no_na_bristol <- GMbC_2500_long_no_na %>% drop_na(bristol_stool_scale)
GMbC_2500_long_no_na_bristol <- GMbC_2500_long_no_na_bristol %>%
        mutate(bristol_stool_scale = as.numeric(bristol_stool_scale))
# Calculate number of non-missing observations per species
observation_counts <- GMbC_2500_long_no_na_bristol %>%
        group_by(species) %>%
        summarise(
                n_obs = sum(!is.na(PTR_value) & !is.na(bristol_stool_scale))
        )

# Filter species with sufficient observations (adjust threshold as needed)
sufficient_data_species <- observation_counts %>%
        filter(n_obs >= 10) %>%  # Example threshold of 10 observations per species
        pull(species)

# Filter GMbC_2500_long_no_na_bristol to retain only species with sufficient data
GMbC_2500_long_filtered <- GMbC_2500_long_no_na_bristol %>%
        filter(species %in% sufficient_data_species)
#perform corrrelation analyis
correlation_result <- GMbC_2500_long_filtered %>%
        group_by(species) %>%
        summarise(
                correlation = cor.test(PTR_value, bristol_stool_scale, method = "spearman", exact = FALSE)$estimate,
                p_value = cor.test(PTR_value, bristol_stool_scale, method = "spearman", exact = FALSE)$p.value
        ) %>%
        mutate(p_adj = p.adjust(p_value, method = "fdr"))

# View the results
print(correlation_result)

# Filter for significant correlations
significant_correlations <- correlation_result %>%
        filter(p_adj < 0.05)
#Only top 30 taxa 
significant_results_top30<-significant_correlations%>%
        inner_join(top_30_taxa,by="species")
#Join with full taxa
full_taxa<-inner_join(significant_correlations,merged_genome,by="species")
#Only keep the first occurances
significant_correlations_taxa <- full_taxa %>%
        distinct(species, .keep_all = TRUE) %>%
        arrange(desc(correlation))
unique(significant_correlations_taxa$phylum)
#Count phylums for results
phylum_counts <- significant_correlations_taxa %>%
        group_by(phylum) %>%
        summarise(count = n()) %>%
        arrange(desc(count))
#save for tsv
write.table(significant_correlations_taxa, file = "significant_correlations_brsitol_taxa_GMbC.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#Visulization
ggplot(significant_correlations, aes(x = reorder(species, correlation), y = correlation, fill = species)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_minimal() +
        labs(
                title = "Spearman Correlations: PTR value vs Bristol stool scale",
                x = "Species",
                y = "Correlation Coefficient",
                fill = "Phylum"
        ) +
        theme(axis.text.y = element_text(face = "italic"))+
        theme(legend.position = "none")

# Save the plot
ggsave("PTR_Bristol_Correlation_Plot_GMbC.jpeg", width = 12, height = 20, dpi = 300)


# Save the plot
ggsave("PTR_Bristol_Correlation_Plot.png", width = 12, height = 10, dpi = 300)
