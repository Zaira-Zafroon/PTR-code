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
KIN_2500 <- read.csv("/Users/bioinfo/Documents/Masters/KIN_2500_out.csv",
                     sep = ",", 
                     header = TRUE, 
                     comment.char = "",
                     check.names = FALSE, 
                     na.strings = "NA")
#head(KIN_2500)

#colnames(KIN_2500)

#Obtain the sepcies name from merged genome
subset_genome <- subset(merged_genome, select = c(user_genome, species))
#Change the column names for merging
names(KIN_2500)[names(KIN_2500) == "log2(PTR):genome_id/sample_id"] <- "genome_id"
names(subset_genome)[names(subset_genome) == "user_genome"] <- "genome_id"
#Clean the sample names
names(KIN_2500) <- gsub("_R1_clean.fastq", "", names(KIN_2500))
KIN_2500 <- KIN_2500 %>%
        left_join(subset_genome, by = "genome_id")
#unique(subset_genome$species)#4335
#Remove the duplicates
KIN_2500 <- KIN_2500 %>%
        distinct()
#Remove the Genome ID column
KIN_2500 <- KIN_2500 %>%
        select(-genome_id)
#Pivot Longer
KIN_2500_long <- KIN_2500 %>%
        pivot_longer(cols = -species, names_to = "SampleID", values_to = "PTR_value")
length(unique(KIN_2500_long$species))#1942 species
#Remove NA 
KIN_2500_long_no_na <- na.omit(KIN_2500_long, cols = "PTR_value")
length(unique(KIN_2500_long_no_na$species))
####relevant taxa####
# CD enriched vector
CD_enriched <- c(
        "s__Enterocloster_clostridioformis",
        "s__Enterocloster_bolteae",
        "s__Enterocloster_clostridioformis_A",
        "s__Enterocloster_aldenensis",
        "s__Escherichia_coli"
)

# CD depleted vector
CD_depleted <- c(
        "s__Faecalibacterium_sp900539945",
        "s__Faecalibacterium_prausnitzii_D",
        "s__Faecalibacterium_longum",
        "s__Faecalibacterium_prausnitzii",
        "s__Ruminococcus_D_bicirculans"
)

# UC enriched vector
UC_enriched <- c(
        "s__Acutalibacter_ornithocaccae",
        "s__Lachnospira_eligens_A",
        "s__Lachnospira_eligens",
        "s__Lachnospira_sp003451515",
        "s__UMGS1071_sp900542375"
)

# UC depleted vector
UC_depleted <- c(
        "s__Bacteroides_caccae",
        "s__Bacteroides_sp905197435",
        "s__Bilophila_wadsworthia",
        "s__Bacteroides_ovatus",
        "s__Bacteroides_xylanisolvens"
)

# IBD enriched vector
IBD_enriched <- c(
        "s__Ruminococcus_B_gnavus",
        "s__Flavonifractor_plautii",
        "s__Bacteroides_fragilis",
        "s__Clostridium_Q_symbiosum",
        "s__Bacteroides_fragilis_A"
)

# IBD depleted vector
IBD_depleted <- c(
        "s__Gemmiger_formicilis",
        "s__ER4_sp000765235",
        "s__Gemmiger_qucibialis",
        "s__Alistipes_A_ihumii",
        "s__Alistipes_A_indistinctus"
)
#join all species
all_species_relevant <- c(
        CD_enriched,
        CD_depleted,
        UC_enriched,
        UC_depleted,
        IBD_enriched,
        IBD_depleted
)
#Create my data frame
KINDRED_relevant_data <- data.frame(
        species = c(CD_enriched, CD_depleted, UC_enriched, UC_depleted, IBD_enriched, IBD_depleted),
        species_status = c(rep("CD_enriched", length(CD_enriched)),
                           rep("CD_depleted", length(CD_depleted)),
                           rep("UC_enriched", length(UC_enriched)),
                           rep("UC_depleted", length(UC_depleted)),
                           rep("IBD_enriched", length(IBD_enriched)),
                           rep("IBD_depleted", length(IBD_depleted)))
)
#Import the metadata file
KIN_metadata <- read.table("metadata_reformat_mgx.tsv", sep = "\t", header = T, 
                           comment.char = "", check.names = F, na.strings = "NA")
KIN_metadata <- KIN_metadata %>%
        rename(SampleID = sample)
head(KIN_metadata)
#Create a merged PTR and metdata file
KIN_2500_long_no_na <- merge(KIN_2500_long_no_na, KIN_metadata, by = "SampleID")
colnames(KIN_2500_long_no_na)
####Bristol stool scale and relationship####
# Drop rows with missing Bristol Stool Scale values
unique(KIN_2500_long_no_na$t13408_Bristol_stool_scale)
# Replace 7777777 with NA in the t13408_Bristol_stool_scale column using subset method
KIN_2500_long_no_na$t13408_Bristol_stool_scale[KIN_2500_long_no_na$t13408_Bristol_stool_scale == 7777777] <- NA
KIN_2500_long_no_na_bristol <- KIN_2500_long_no_na %>% drop_na(t13408_Bristol_stool_scale)
KIN_2500_long_no_na_bristol <- KIN_2500_long_no_na_bristol %>%
        mutate(bristol_stool_scale = as.numeric(t13408_Bristol_stool_scale))
unique(KIN_2500_long_no_na_bristol$bristol_stool_scale)
# Calculate number of non-missing observations per species
observation_counts <- KIN_2500_long_no_na_bristol %>%
        group_by(species) %>%
        summarise(
                n_obs = sum(!is.na(PTR_value) & !is.na(bristol_stool_scale))
        )

# Filter species with sufficient observations (adjust threshold as needed)
sufficient_data_species <- observation_counts %>%
        filter(n_obs >= 10) %>%  # Example threshold of 10 observations per species
        pull(species)

# Filter KIN_2500_long_no_na_bristol to retain only species with sufficient data
KIN_2500_long_filtered <- KIN_2500_long_no_na_bristol %>%
        filter(species %in% sufficient_data_species)
#perform corrrelation analyis
correlation_result <- KIN_2500_long_filtered %>%
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
        filter(p_adj < 0.05)#no significant correlation
significant_correlations_p_value <- correlation_result %>%
        filter(p_value < 0.05)
unique(significant_correlations_p_value$species)

#Visulization
ggplot(significant_correlations_p_value, aes(x = reorder(species, correlation), y = correlation, fill = species)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_minimal() +
        labs(
                title = "Spearman Correlations: PTR value vs Bristol stool scale",
                x = "Species",
                y = "Correlation Coefficient",
                fill = "Species"
        ) +
        theme(axis.text.y = element_text(face = "italic")) +
        theme(legend.position = "none")

# Save the plot
ggsave("PTR_Bristol_Correlation_Plot_KIN.jpeg", width = 10, height = 6, dpi = 300)


