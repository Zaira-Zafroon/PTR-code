#####library####
library(dplyr)
library(tidyverse)
library(purrr)
library(rstatix)
library(broom)
#####Import files####
#Import the previously coded merged genome for isolates and MAG assignment
merged_genome<-read.table("merged_genome.tsv",
                          sep = "\t", 
                          header = T,  comment.char = "",
                          check.names = F, na.strings = "NA")
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
unique(GMbC_2500_long_no_na$lifestyle)
GMbC_2500_long_no_na %>%
        group_by(species) %>%
        summarize(unique_levels = n_distinct(lifestyle)) %>%
        filter(unique_levels == 1)#846 

#Only keep data where there are two levels of lifestyle
filtered_data <- GMbC_2500_long_no_na %>%
        group_by(species) %>%
        filter(n_distinct(lifestyle) > 1) %>%
        ungroup()
#hist(filtered_data$PTR_value, main = "Histogram of PTR Values", xlab = "PTR Values")#not normally ditibuted

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
significant_species_associaition<-significant_association_results_lifestyle %>%
        select(species,p.adj)
filtered_data_sig <- filtered_data %>%
        inner_join(significant_species_associaition, by= "species") 
#Find the effect sizes for PTR values and lifestyle
significant_effect_size_lifestyle<-filtered_data_sig %>%
        nest(data= -species)%>%#creates a nested data frame where the species column not netsed,rest nested
        #mutate will create new column
        #map function allows you to take a function/test and reiterate it over data
        mutate(test=map(.x=data, ~wilcox_effsize(PTR_value~lifestyle,data=.x)))%>%
        #wilcox effect test performed on specified columns inside the tibble per species
        #the output is not a htest but a direct table and can be unnested
        unnest(test)#makes the test tiblble directly part of data frame
colnames(significant_effect_size_lifestyle)
#Merge the dataset with p adjusted values        
significant_species<-merge(significant_species_associaition,significant_effect_size_lifestyle,by="species")
colnames(significant_species)
#For tsv files remove data column
significant_species_table<-significant_species %>%
        select(-data)
colnames(significant_species_table)
#Write the file
write.table(significant_species_table, file = "significant_species_table.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Filter for large effect sizes 
large_effect_size_species<-significant_species %>%
        filter(magnitude=="large")

#Do a clean and prep it for a visual
large_effect_size_association<-filtered_data %>%
        inner_join(large_effect_size_species, by = "species") %>%
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
        inner_join(large_effect_size_species, by = "species") %>%
        ggplot(aes(x = PTR_value, y = lifestyle, color = country, shape = urbanism)) +
        geom_jitter(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3),
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

ggsave(filename = "large_effect_size_association_higher_res.png", plot = large_effect_size_association_higher_res, width = 12, height = 10, units = "in", dpi = 300)

# If you want to adjust the order of the facets
# large_effect_size_association + facet_wrap(~ fct_reorder(species, PTR_value, .fun = median), scales = "free_y", ncol = 2)
#Find the moderate effect size species
moderate_effect_size_species<-significant_species %>%
        filter(magnitude=="moderate")
#Visualization

moderate_effect_size_species_lifestyle<-filtered_data %>%
        inner_join(moderate_effect_size_species, by="species") %>%
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
        filename = "moderate_effect_size_species_lifestyle.png", 
        plot = moderate_effect_size_species_lifestyle, 
        device = "png", 
        width = 15, 
        height = 28, 
        units = "in", 
        dpi = 300
)