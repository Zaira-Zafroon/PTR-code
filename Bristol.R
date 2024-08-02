####data import and cleaning####

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