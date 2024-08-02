####library load ####
library(FSA)

library(ggnewscale)

library(partR2)

library(ggplot2)

library(tidyverse)

library(lme4)

library(ggplot2)

library(dplyr)

library(mia)

library(purrr)

library(tibble)

library(ggplot2)

library(gganimate)

library(transformr)

library(plotly)

library(tidyr)

library(stringr)

library(broom)

library(purrr)

library(rstatix)

library(ggstatsplot)

library(ggtext)

library(lme4)

###Import files ####
#import taxonomy files
#Import the SGB representative as genomes are named by clusters
GloHugg_rep<- read.table("GloHuGG.GTDBr214.cluster_final_tax.tsv", sep = "\t", 
                         header = T,  comment.char = "",
                         check.names = F, na.strings = "NA")
head(GloHugg_rep$classification)


#Separate the data into columns according to taxonomy

GloHugg_rep <- separate(GloHugg_rep, classification, 
                        into = c("domain", "phylum", "class", "order", "family", "genus", "species"), 
                        sep = ";", extra = "merge")
unique(GloHugg_rep$species)


#Modify the taxon names in the GloHugg table to match
GloHugg_rep$species<- gsub(" ", "_",GloHugg_rep$species)

#The genome information-genome source 
GloHugg_genome<-read.table("/Users/bioinfo/Documents/Masters/GloHuGG.all_genomes.tsv",
                           sep = "\t", 
                           header = T,  comment.char = "",
                           check.names = F, na.strings = "NA")
colnames(GloHugg_genome)
GloHugg_genome<- GloHugg_genome %>% 
        select(-Domain, -Score, -Completeness, -Contamination,
               -cluster_95_final,-cluster_97_final,- cluster_99_final,
               -rep_preclust)
GloHugg_genome<- rename(GloHugg_genome, user_genome = final)
####Creating Merged Database####
merged_genome<- merge(GloHugg_rep, GloHugg_genome, by = "user_genome")
dim(merged_genome)
unique(merged_genome$species)
write.table(merged_genome, file = "merged_genome.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
####Import Merged genome instead of running the code####
merged_genome<-read.table("merged_genome.tsv",
                           sep = "\t", 
                           header = T,  comment.char = "",
                           check.names = F, na.strings = "NA")

#Add new column specifying the source
library(dplyr)

# Assuming your dataset is called 'df'
merged_genome <- merged_genome %>%
        mutate(source = case_when(
                grepl("GMBC|^BN10", genome) ~ "isolates",
                !grepl("GMBC|^BN10", genome) ~ "assembled genomes"
        ))

k<-unique(merged_genome$species)#4335
####Obtaining Unique Species List####
#Find unique species that are sourced from isolates
unique_species_isolates <- merged_genome %>%
        group_by(species) %>%
        filter(source == "isolates") %>%
        slice(1) %>%
        ungroup() #420 species from isolates 
#Find unique species that are sourced from assmeblies
unique_species_assembled_genome <- merged_genome %>%
        group_by(species) %>%
        filter(source == "assembled genomes") %>%
        slice(1) %>%
        ungroup()    #4303 from assembled genome     
#Merge them to find a list of unique speices
unique_species <- unique_species_isolates %>% #kept as is
        bind_rows(
                unique_species_assembled_genome %>% #takes this and create a new data frame that is non reducndant
                        anti_join(unique_species_isolates, by = "species")
        )

k<-unique(unique_species$species)
length(k)
#Clean the data frame 
unique_species <- unique_species %>% 
        select(-domain, -phylum, -class, -order, -family, -genus,-genome,-source)

####For copy Script bash####
files_with_extension <- paste0(unique_species$user_genome, ".fna")
files_with_extension_source <- paste0(" /work_beegfs/ikmb_repository/shared/microbiome/processed_datasets/collections/GloHuGG/SGB_representatives/genomes/",files_with_extension)
files_with_extension_destination <- paste0(files_with_extension_source," /work_beegfs/sukmb625/coptr/GMbC/ref-db/GMbC_ref")
write(files_with_extension_destination, file = "coptr_ref.txt", sep = "\n")

####ACTUAL STATS START HERE####
#####Import files####
#Import GMbC_metadata
GMbC_metadata<-read.table("metadata_gmbc_bn10_complete.tsv", sep = "\t", header = T, 
                          comment.char = "",check.names = F, na.strings = "NA")
head(GMbC_metadata)
#Create Sample Column
GMbC_metadata<- tibble::rownames_to_column(GMbC_metadata, "sample")

unique(GMbC_metadata$bristol_stool_scale)
#Import COPTR results
GMbC_2500<-read.csv("/Users/bioinfo/Documents/Masters/GMbC_2500_out.csv",
                      sep = ",", 
                      header = TRUE, 
                      comment.char = "",
                      check.names = FALSE, 
                      na.strings = "NA")
head(GMbC_2500)

colnames(GMbC_2500)
unique(GMbC_metadata$iga)
#Obtain the sepcies name from merged gebome
subset_genome <- subset(merged_genome, select = c(user_genome,species))
#Change the column names for merging
names(GMbC_2500)[names(GMbC_2500) == "log2(PTR):genome_id/sample_id"] <- "genome_id"
names(subset_genome)[names(subset_genome) == "user_genome"] <- "genome_id"
#Clean the sample names
names(GMbC_2500) <- gsub("_R1_clean.fastq", "", names(GMbC_2500))
GMbC_2500 <- GMbC_2500 %>%
        left_join(subset_genome, by = "genome_id")
unique(subset_genome$species)
#Remove the duplicates
GMbC_2500 <- GMbC_2500 %>%
        distinct()
#Remove the Genome ID column
GMbC_2500 <- GMbC_2500 %>%
        select(-genome_id)
#Pivot Longer
GMbC_2500_long <- GMbC_2500 %>%
        pivot_longer(cols = -species, names_to = "sample", values_to = "PTR_value")
length(unique(GMbC_2500_long$species))#1942 species
#Remove NA 
GMbC_2500_long_no_na <- na.omit(GMbC_2500_long, cols = "PTR_value")
length(unique(GMbC_2500_long_no_na$species))#1942 species
#Merge the metadata with PTR values in the long format
GMbC_2500_long_no_na<-merge(GMbC_2500_long_no_na,GMbC_metadata,by="sample")
#Hitogram of relationship
hist(GMbC_2500_long_no_na$PTR_value, main = "Histogram of PTR Values", xlab = "PTR Values")
qqnorm(GMbC_2500_long_no_na$PTR_value)#Not normally distributed


####Basic Studies on mean PTR per individual####

summary_stats <- GMbC_2500_long_no_na %>%
        group_by(species) %>%
        summarise(
                Mean = mean(PTR_value, na.rm = TRUE),
                Median = median(PTR_value, na.rm = TRUE),
                SD = sd(PTR_value, na.rm = TRUE),
                Min = min(PTR_value, na.rm = TRUE),
                Max = max(PTR_value, na.rm = TRUE),
                Count = n()
        )

#Assciation study results per species PTR with lifestyle 


#Check data distribution for values
hist(GMbC_2500_long$PTR_value, main = "Histogram of PTR Values", xlab = "PTR Values")#not normally ditibuted

####Do the species relative abundance followup toc ehck PTRS are uncorelated with RA####
????
###lifestyle#####
#Remove NA values from the data 
GMbC_2500_long_no_na <- na.omit(GMbC_2500_long, cols = "PTR_value")
unique(GMbC_2500_long_no_na$lifestyle)
GMbC_2500_long_no_na %>%
        group_by(species) %>%
        summarize(unique_levels = n_distinct(lifestyle)) %>%
        filter(unique_levels == 1)#879 

#Only keep data where there are two levels of lifestyle
filtered_data <- GMbC_2500_long_no_na %>%
        group_by(species) %>%
        filter(n_distinct(lifestyle) > 1) %>%
        ungroup()
hist(filtered_data$PTR_value, main = "Histogram of PTR Values", xlab = "PTR Values")#not normally ditibuted



####Overview Data ####
#GMbC_2500_long the long format data
#the ranked taxa is from the RefGenomeGMbC_new.R script
ranked_taxa<- read.table("ranked_taxa.tsv", sep = "\t", 
                         header = T,  comment.char = "",
                         check.names = F, na.strings = "NA")
top_30_taxa <- ranked_taxa %>% 
        slice(1:30)
# Renaming the column in the data frame
top_30_taxa  <- top_30_taxa  %>%
        rename(species = taxon)
observations_per_species <- GMbC_2500_long_no_na %>%
        group_by(species) %>%
        summarise(num_observations = n())
median_species <- median(observations_per_species$num_observations)
summary_stats <- GMbC_2500_long_no_na %>%
        inner_join(top_30_taxa, by = "species") %>%
        group_by(species) %>%
        summarise(
                Mean = mean(PTR_value, na.rm = TRUE),
                Median = median(PTR_value, na.rm = TRUE),
                SD = sd(PTR_value, na.rm = TRUE),
                Min = min(PTR_value, na.rm = TRUE),
                Max = max(PTR_value, na.rm = TRUE),
                Count = n()
        )
print(summary_stats)
mean_PTR_value_distribution <- ggplot(summary_stats, aes(x = species, y = Mean, fill = species)) +
        geom_point(shape = 21, size = 3) +  # Use shape 21 for filled circles
        labs(title = "Distribution of Mean PTR Values by Taxon",
             x = "Taxon", y = "Mean PTR Value") +
        theme_classic() +  # Apply classic theme
        theme(axis.text.x = element_text(angle = 65, hjust = 1))
        
# Save the plot as a JPEG with high resolution
ggsave(
        filename = "mean_PTR_value_distribution.jpeg", 
        plot = mean_PTR_value_distribution, 
        device = "jpeg", 
        width = 10, 
        height = 8, 
        units = "in", 
        dpi = 300
)

GMbC_2500_long_no_na %>%
        inner_join(top_30_taxa, by = "species") %>%
        mutate(
                taxonomy = str_replace(species, "(.*)", "*\\1*"),
                taxonomy = str_replace(taxonomy, "\\*(.*)_unclassified\\*", "Unclassified<br>*\\1*")
        ) %>%
        ggplot(aes(x = PTR_value, y = taxonomy, fill = species)) +
        geom_boxplot(color = "black", outlier.color = "black") +
        labs(x = "PTR_value", y = NULL)  +
        theme_classic() +
        theme(axis.text.y = element_markdown(),
              legend.position = "none")  # Remove legend
PTR_value_distribution<-GMbC_2500_long_no_na %>%
        inner_join(top_30_taxa, by = "species") %>%
        mutate(
                taxonomy = str_replace(species, "(.*)", "*\\1*"),
                taxonomy = str_replace(taxonomy, "\\*(.*)_unclassified\\*", "Unclassified<br>*\\1*"),
                lifestyle = factor(lifestyle, levels = c("non_industrialized", "industrialized"))
        ) %>%
        ggplot(aes(x = PTR_value, y = taxonomy, fill = lifestyle)) +
        geom_boxplot(color = "black", outlier.color = "black") +
        labs(x = "PTR_value", y = NULL)  +
        theme_classic() +
        theme(axis.text.y = element_markdown(),
              legend.position = "none")
# Save the plot as a JPEG with high resolution
ggsave(
        filename = "PTR_value_distribution.jpeg", 
        plot = PTR_value_distribution, 
        device = "jpeg", 
        width = 10, 
        height = 8, 
        units = "in", 
        dpi = 300
)

#Find the significant association results with lifestlye and PTR
significant_association_results_lifestyle<-filtered_data %>%
        nest(data= -species)%>%#creates a nested data frame where the species column not netsed,rest nested
#mutate will create new column
#map function allows you to take a function/test and reiterate it over data
        mutate(test=map(.x=data, ~wilcox.test(PTR_value~lifestyle,data=.x)%>%tidy))%>%
#wilcox test performed on specified columns inside the tibble per species
#test results are a htest object that need to be tidied
#Using broom package for that and the function tidy to make it a tibble
#Tidy has to be written in the map part as written now
unnest(test)%>%#makes the test tiblble directly part of data frame
mutate(p.adj=p.adjust(p.value,method="fdr"))%>%
        filter(p.adj<0.05)
#Streamline the columns
significant_genera_associaition<-significant_association_results_lifestyle %>%
        select(species,p.adj)

#Find the effect sizes for PTR values and lifestyle
significant_effect_size_lifestyle<-filtered_data %>%
        nest(data= -species)%>%#creates a nested data frame where the species column not netsed,rest nested
        #mutate will create new column
        #map function allows you to take a function/test and reiterate it over data
        mutate(test=map(.x=data, ~wilcox_effsize(PTR_value~lifestyle,data=.x)))%>%
        #wilcox effect test performed on specified columns inside the tibble per species
        #the output is not a htest but a direct table and can be unnested
        unnest(test)#makes the test tiblble directly part of data frame

#Merge the dataset with p adjusted values        
significant_genera<-merge(significant_genera_associaition,significant_effect_size_lifestyle,by="species")
large_effect_size_genera<-significant_genera %>%
        filter(magnitude=="large")

#Do a clean and prep it for a visual
large_effect_size_association<-filtered_data %>%
        inner_join(large_effect_size_genera, by = "species") %>%
        ggplot(aes(x = PTR_value, y = species, colour = lifestyle, fill = lifestyle)) +
        geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3),
                    shape = 21) + #cosmetics
        #colour and fill stating both of them together with shape 21 allows for a inner circle with a boudnary line
        stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
                     geom = "pointrange",
                     position = position_dodge(width = 0.8),
                     color = "black", show.legend = FALSE)#that shape commanda before allows us to make it black
#Clear differences in the median values of the species
# Save the plot as a JPEG with high resolution
ggsave(
        filename = "large_effect_size_association_with_lifestyle.png", 
        plot = large_effect_size_association, 
        device = "png", 
        width = 10, 
        height = 8, 
        units = "in", 
        dpi = 300
)
#Higher Resoultion Visualization

# Define a 12-color palette
color_palette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
                   "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
                   "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
large_effect_size_association_higher_res <- filtered_data %>%
        inner_join(large_effect_size_genera, by = "species") %>%
        ggplot(aes(x = PTR_value, y = lifestyle, color = country, shape = urbanism)) +
        geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3),
                    size = 3) +
        stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
                     geom = "pointrange",
                     position = position_dodge(width = 0.8),
                     color = "black", show.legend = FALSE) +
        facet_wrap(~ species, scales = "free_y", ncol = 2) +
        theme_classic() +
        theme(
                strip.background = element_rect(fill = "lightgrey"),
                strip.text = element_text(face = "italic"),
                legend.position = "bottom",
                axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(x = "PTR Value", y = "Lifestyle", color = "Country", shape = "Urbanism") +
        scale_color_manual(values = color_palette) +
        scale_shape_manual(values = c(16, 17)) 
ggsave(
        filename = "large_effect_size_association_with_lifestyle_country_urbanism.png", 
        plot = large_effect_size_association_higher_res, 
        device = "png", 
        width = 28, 
        height = 12, 
        units = "in", 
        dpi = 300
)
# If you want to adjust the order of the facets
# large_effect_size_association + facet_wrap(~ fct_reorder(species, PTR_value, .fun = median), scales = "free_y", ncol = 2)
#Find the moderate effect size genera
moderate_effect_size_genera<-significant_genera %>%
        filter(magnitude=="moderate")
#Visualization

moderate_effect_size_genera_lifestyle<-filtered_data %>%
        inner_join(moderate_effect_size_genera, by="species") %>%
        mutate(
                taxonomy = str_replace(species, "(.*)", "*\\1*"),#to make the species italic first put stars at the side of the names 
                taxonomy = str_replace(taxonomy, "\\*(.*)_unclassified\\*", "Unclassified<br>*\\1*"),#to make the species italic when its unclassified
                #the upper two lines put the stars on the spcies name
                #when i apply ggtext it would do its job
               lifestyle = factor(lifestyle, levels = c("non_industrialized", "industrialized"))#order the way I want
        ) %>%
        ggplot(aes(x = PTR_value, y = taxonomy, color = lifestyle, fill = lifestyle)) +
        geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5), shape = 21) +
        stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), geom = "pointrange", 
                     position = position_dodge(width = 0.8), color = "black", show.legend = FALSE) +
       scale_color_manual(NULL, #for the colour
                           breaks = c("industrialized", "non_industrialized"), 
                           values = c("gray", "dodgerblue"), 
                           labels = c("industrialized", "non_industrialized")) +
        scale_fill_manual(NULL, #for the fill
                          breaks = c("industrialized", "non_industrialized"), 
                          values = c("gray", "dodgerblue"), 
                          labels = c("industrialized", "non_industrialized")) +
        labs(x = "PTR_value", y = NULL) + #the totle
        theme_classic() + #no grid line 
        theme(axis.text.y = element_markdown())# themarkdown for italicized text in y axis
ggsave(
        filename = "moderate_effect_size_genera_lifestyle.png", 
        plot = moderate_effect_size_genera_lifestyle, 
        device = "png", 
        width = 10, 
        height = 8, 
        units = "in", 
        dpi = 300
)
###Bristol Stool Scale###
#Take the without missing PTR values columns
"#####Bristol score scale#####
filtered_bristo_check<-GMbC_2500_long_no_na %>%
        group_by(species) %>%
        summarize(unique_levels = n_distinct(bristol_stool_scale)) %>%
        filter(unique_levels == 1)#22
filtered_bristol <- GMbC_2500_long_no_na %>%
        group_by(species) %>%
        filter(n_distinct(bristol_stool_scale) > 1) %>%
        ungroup()
filtered_bristol <- as.data.frame(filtered_bristol)
bristol_signifiance<-filtered_bristol %>% #more than one group pf bristol score per species
        nest(data=-species)%>%#cretae a second column as a list o table
        mutate(experiment_tests=map(.x=data,
                                    ~kruskal.test(PTR_value~bristol_stool_scale,data=.x)%>% 
                                    tidy())) %>% 
# tidy takes htest and make it into data frame/tiblle and it only works with base R kruskal test
                                unnest(experiment_tests)%>% #get the results out of dataframe
                                mutate(p.adj=p.adjust(p.value,method="BH"))%>% #correction for multiple test
                                 filter(p.adj<0.05)%>% #filter data
        select(species,data,p.adj)#data fields
#checked if each species had a signfcant difference in PTR between all seven bristol score groups 

#Kept the nested data frame so that we can do the pairwise test directly

significant_bristol_pairwise<-bristol_signifiance %>%
        mutate(pairwise_tests = map(.x=data,
                                    ~pairwise.wilcox.test(x=.x$PTR_value,
                                                          g=.x$bristol_stool_scale,
                                                          p.adjust.method = "BH") %>%
                                            tidy())) %>%
        unnest(pairwise_tests) %>%
        filter(p.value < 0.05) %>%
        select(species,data, group1, group2, p.value) %>%
        ungroup()                 
unique(significant_bristol_pairwise$species)
#overview data
bristol_species<-unique(significant_bristol_pairwise$species)#33
#Visulization
bristol_overview <- filtered_data %>%
        filter(species %in% bristol_species) %>%
        mutate(bristol_stool_scale = factor(bristol_stool_scale)) %>% # Convert to factor
        group_by(species) %>%
        ggplot(aes(x = PTR_value, y = species, color = bristol_stool_scale, fill = bristol_stool_scale)) +
        geom_boxplot() +
        stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), geom = "pointrange",
                     position = position_dodge(width = 0.8), color = "black", show.legend = FALSE) +
        scale_color_manual(NULL,
                           breaks = c("1", "2", "3", "4", "5", "6", "7"),
                           values = c("#808080", "#1E90FF", "#FF0000", "#FFA500", "#008000", "#FFC0CB", "#FFFF00"),
                           labels = c("1", "2", "3", "4", "5", "6", "7")) +
        scale_fill_manual(NULL,
                          breaks = c("1", "2", "3", "4", "5", "6", "7"),
                          values = c("#808080", "#1E90FF", "#FF0000", "#FFA500", "#008000", "#FFC0CB", "#FFFF00"),
                          labels = c("1", "2", "3", "4", "5", "6", "7")) +
        labs(x = "PTR_value", y = NULL) +
        theme_classic() +
        theme(axis.text.y = element_markdown())
ggsave(
        filename = "basic_bristol_plot.png", 
        plot = bristol_overview, 
        device = "png", 
        width = 10, 
        height = 8, 
        units = "in", 
        dpi = 300
)

plots_list <- bristol_signifiance %>%
        mutate(plot = map(data, ~ ggbetweenstats(
                data = .x,
                x = bristol_stool_scale,
                y = PTR_value,
                type = "nonparametric",
                plot.type = "boxviolin",
                pairwise.comparisons = TRUE,
                pairwise.display = "significant",
                p.adjust.method = "bonferroni",
                pairwise.tests = "wilcox.test"
        )))
species_plots <- plots_list %>%
        select(species, plot)""

####iga levels and calpro levels####
GMbC_2500_long_no_na
#Visualize the IgA Levels
# Spearman correlation calculation for each species
GMbC_2500_long_no_na_iga <- GMbC_2500_long_no_na %>% drop_na(iga)
# Calculate number of non-missing observations per species
observation_counts <- GMbC_2500_long_no_na_iga %>%
        group_by(species) %>%
        summarise(
                n_obs = sum(!is.na(PTR_value) & !is.na(iga))
        )

# Filter species with sufficient observations (adjust threshold as needed)
sufficient_data_species <- observation_counts %>%
        filter(n_obs >= 10) %>%  # Example threshold of 10 observations per species
        pull(species)

# Filter GMbC_2500_long_no_na_iga to retain only species with sufficient data
GMbC_2500_long_filtered <- GMbC_2500_long_no_na_iga %>%
        filter(species %in% sufficient_data_species)

# Perform Spearman correlation tests for filtered data
correlation_results_iga <- GMbC_2500_long_filtered %>%
        group_by(species) %>%
        summarise(
                correlation = cor(PTR_value, iga, method = "spearman", use = "pairwise.complete.obs"),
                p.value = cor.test(PTR_value, iga, method = "spearman", exact = FALSE, na.action = na.omit)$p.value
        ) %>%
        mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
        filter(p.adj < 0.05)#34 species
# Visualization for IgA Levels vs. PTR Values

GMbC_2500_long_filtered %>%
        inner_join(correlation_results_iga, by = "species") %>%
        inner_join(top_30_taxa, by = "species") %>%
        group_by(species) %>%
        ggplot(aes(x = iga, y = PTR_value)) +
        geom_point(aes(color = urbanism)) +  # Add color aesthetic here
        geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~species, scales = "free") +
        labs(x = "IgA Levels", y = "PTR Values", color = "Country") +  # Add legend title
        theme_bw() +
        theme(legend.position = "bottom") 

        
####Visulization of Calpro levels####
GMbC_2500_long_no_na_calpro <- GMbC_2500_long_no_na %>% drop_na(calpro)

# Calculate number of non-missing observations per species
observation_counts <- GMbC_2500_long_no_na_calpro %>%
        group_by(species) %>%
        summarise(
                n_obs = sum(!is.na(PTR_value) & !is.na(calpro))
        )

# Filter species with sufficient observations (adjust threshold as needed)
sufficient_data_species <- observation_counts %>%
        filter(n_obs >= 10) %>%  # Example threshold of 10 observations per species
        pull(species)

# Filter GMbC_2500_long_no_na_calpro to retain only species with sufficient data
GMbC_2500_long_filtered <- GMbC_2500_long_no_na_calpro %>%
        filter(species %in% sufficient_data_species)

correlation_results_calpro <- GMbC_2500_long_filtered %>%
        group_by(species) %>%
        summarise(
                correlation = cor(PTR_value, calpro, method = "spearman", use = "pairwise.complete.obs"),
                p.value = cor.test(PTR_value, calpro, method = "spearman", exact = FALSE, na.action = na.omit)$p.value
        ) %>%
        mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
        filter(p.adj < 0.05)

GMbC_2500_long_filtered %>%
        inner_join(correlation_results_calpro, by = "species") %>%
        group_by(species) %>%
        ggplot(aes(x = iga, y = PTR_value)) +
        geom_point(aes(color = urbanism)) +  # Add color aesthetic here
        geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~species, scales = "free") +
        labs(x = "Calprotectin Levels", y = "PTR Values", color = "Country") +  # Add legend title
        theme_bw() +
        theme(legend.position = "bottom") 

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

# Perform linear regression for filtered data
regression_results <- GMbC_2500_long_filtered %>%
        group_by(species) %>%
        nest() %>%
        mutate(
                model = map(data, ~ lm(PTR_value ~ bristol_stool_scale, data = .)),
                tidied = map(model, tidy),
                unnested = map(tidied, ~ filter(.x, term == "bristol_stool_scale")),
                results = map(unnested, ~ .x %>%
                                      select(term, estimate, std.error, statistic, p.value))  # Adjust select() here
        ) %>%
        unnest(results) %>%
mutate(p.adj = p.adjust(p.value, method = "fdr"))%>%
        filter(p.adj<0.05)%>%
        inner_join(top_30_taxa,by="species")
GMbC_2500_long_filtered %>%
        inner_join(regression_results, by = "species") %>%
        ggplot(aes(x = bristol_stool_scale, y = PTR_value)) +
        geom_point(aes(color = urbanism), alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "blue") +
        facet_wrap(~species, scales = "free") +
        labs(x = "Bristol Stool Scale", y = "PTR Values", color = "Urbanism") +
        theme_bw() +
        theme(legend.position = "bottom")
###RA and PTR data####
GMbC_abundance_table<- read.table("GMbC.GTDBr214.species.tpm.out", sep = "\t", 
                                 header = T,  comment.char = "",
                                 check.names = F, na.strings = "NA")

unique(GMbC_abundance_table$set)#only GMbC

#remove the set column 
GMbC_abundance_table <- subset(GMbC_abundance_table, select = -set)
#Clean the GMbC_abundance_table using metadata
GMbC_abundance_table <- GMbC_abundance_table[GMbC_abundance_table$sample %in% GMbC_metadata$sample, ]
unique(GMbC_abundance_table$set)#only GMbC

#remove the set column 
GMbC_abundance_table <- subset(GMbC_abundance_table, select = -set)
#Long format the data
GMbC_abundance_long<-GMbC_abundance_table%>%
pivot_longer(cols = -sample, names_to = "species", values_to = "tpm")
#log transform the tpm
GMbC_abundance_long$log2_tpm <- log2(GMbC_abundance_long$tpm + 1)
#Filter data
species_counts <- GMbC_abundance_long %>%
        group_by(species) %>%
        tally() %>%
        filter(n >= 20) %>%
        pull(species)

filtered_abundance <- GMbC_abundance_long %>%
        filter(species %in% species_counts)

filtered_ptr <- GMbC_2500_long_no_na %>%
        filter(species %in% species_counts)

#Combine the data frame
combined_data <- left_join(filtered_ptr, filtered_abundance, by = c("species", "sample"))
species_corrs <- combined_data %>%
        group_by(species) %>%
        filter(sum(!is.na(PTR_value) & !is.na(log2_tpm)) >= 20) %>%
        filter(!is.infinite(PTR_value) & !is.infinite(log2_tpm) & !is.na(PTR_value) & !is.na(log2_tpm)) %>%
        summarise(correlation = cor(PTR_value, log2_tpm, method = "spearman", use = "pairwise.complete.obs"),
                  p.value = cor.test(PTR_value, log2_tpm, method = "spearman", exact = FALSE, na.action = na.omit)$p.value) %>%
        mutate(p.adj = p.adjust(p.value, method = "fdr"))%>%
        filter(p.adj<0.05)
#Attach_degree of relationship
association_data <- species_corrs %>%
        mutate(degree_of_relationship = case_when(
                abs(correlation) < 0.2  ~ "Very weak",
                abs(correlation) >= 0.2 & abs(correlation) < 0.4 ~ "Weak",
                abs(correlation) >= 0.4 & abs(correlation) < 0.6 ~ "Moderate",
                abs(correlation) >= 0.6 & abs(correlation) < 0.8 ~ "Strong",
                abs(correlation) >= 0.8 & abs(correlation) <= 1.0 ~ "Very strong"
        ))


# Define colors for positive and negative correlations
color_palette <- scale_fill_manual(values = c(
        "Very weak" = "lightblue", 
        "Weak" = "blue", 
        "Moderate" = "orange", 
        "Strong" = "red"
), name = "Degree of Relationship")

# Create a plot using ggplot2
heatmap_plot <- ggplot(heatmap_data, aes(x = 1, y = species)) +
        geom_tile(aes(fill = degree_of_relationship), color = "white") +
        geom_point(aes(fill = degree_of_relationship, size = abs(correlation)), shape = 21, color = "black") +
        scale_size_continuous(range = c(2, 5)) +  # Adjust size range as needed
        color_palette +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 7)) +
        labs(x = NULL, y = "Species", title = "Spearman Correlation between PTR and Log-Transformed Relative Abundance")

print(heatmap_plot)
