library(ggpattern)

library(dendextend)

library(ggdendro)

library(cowplot)

# Read OTU table for SC100 screening
ot_scree_filtered <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/scree_ot_filtered.csv",
                        header = T,
                        row.names=1)

# Minmax transform
ot_scree_filtered_norm <- transform_feature_table(ot_scree_filtered, transform_method = "min_max")
# Get clustering order
ordered_samples_cluster <- order_samples_by_clustering(ot_scree_filtered_norm)

# Generate df long format for Species abundance data
ot_scree_filtered2 <- ot_scree_filtered %>% rownames_to_column(var = "Species")
df_otu_long <- ot_scree_filtered2 %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")

# Update sample factor levels in the long-format data for ggplot
df_otu_long$Sample <- factor(df_otu_long$Sample, levels = ordered_samples_cluster)


##### Now lets work with the strain data
strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA31", col_names = TRUE)

# Define strain numbers based on the position of each strain within the species
strain_data2 <- strain_data %>%
  mutate(Strain_Number = rep(1:3, times = 10)) # 1, 2, 3 for each species

# Replace `1`s with the strain number for each sample, keeping `0`s unchanged
for (sample_col in colnames(strain_data2[4:ncol(strain_data2)-1])) {
  strain_data2[[sample_col]] <- ifelse(strain_data2[[sample_col]] == 1, strain_data2$Strain_Number, 0)
}

# Drop the helper "Strain_Number" column to get the desired output
df_final <- strain_data2 %>% select(-Strain_Number)

# Collapse the information by species
df_collapsed <- df_final %>%
  mutate(Species = sapply(strsplit(Species, " "), function(x) paste(x[1:2], collapse = " "))) %>% # Extract species name
  group_by(Species) %>%
  summarise(across(starts_with("SC"), max)) %>% # Take max per sample to represent strain
  ungroup()

new_row <- c("Staphylococcus aureus", rep(1, ncol(df_collapsed) - 1))
new_row_df <- as.data.frame(t(new_row), stringsAsFactors = FALSE)
colnames(new_row_df) <- colnames(df_collapsed)

df_collapsed2 <- rbind(df_collapsed, new_row_df)

##### Generate table for 
df_otu_long2 <- df_collapsed2 %>% # this might be df_strain_long
  pivot_longer(-Species, names_to = "Sample", values_to = "Strain")


### Join the two long dataframes

df_joined <- dplyr::full_join(df_otu_long, df_otu_long2, by = c("Species", "Sample"))

df_joined_filtered <- df_joined %>%
  filter(!is.na(Abundance) & Abundance != 0)

df_joined_filtered <- df_joined_filtered %>%
  filter(!is.na(Strain) & Strain != 0)

df_joined_filtered <- df_joined_filtered %>%
  mutate(Strain = factor(
    case_when(
      Strain == 1 ~ "Strain 1",
      Strain == 2 ~ "Strain 2",
      Strain == 3 ~ "Strain 3"
    )
  ))

# Update sample factor levels in the long-format data for ggplot
df_joined_filtered$Sample <- factor(df_joined_filtered$Sample, levels = ordered_samples_cluster)

####### Now plotting

colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1","brown1", "olivedrab3", "darkorange3", "#23001E","hotpink" )


# this is the one, do not touch
p1 <- ggplot(data = df_joined_filtered, aes(x = Sample, y=Abundance)) +
  geom_bar_pattern(aes(fill = Species, pattern = Strain, pattern_density = Strain),
                   position = "fill",
                   stat="identity",
                   show.legend = TRUE,
                   pattern_color = "white",
                   pattern_fill = "white",
                   pattern_angle = 45,
                   pattern_spacing = 0.025) +
  ggplot2::scale_fill_manual(values=colours_vec) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                 axis.text.y = ggplot2::element_text(size=12),
                 legend.text = element_text(size=12)) +
  guides(pattern = guide_legend(override.aes = list(fill = "black")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  scale_pattern_manual(values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")) +
  scale_pattern_density_manual(values = c(0, 0.2, 0.05))

p1

ggsave(
  "C:/Users/marce/Desktop/barplot_scree.pdf",
  plot = p1,
  device = "pdf",
  width = 20,
  height = 10, scale = 0.5
)

p2 <- ggdendrogram(hc, rotate = 0,
                   leaf_labels = F) +
  theme(plot.margin = margin(t = 40,  # Top margin
                             r = 0,  # Right margin
                             b = -25,  # Bottom margin
                             l = 0)) + # Left margin
  scale_y_continuous(expand = expansion(add = c(1, 1)), labels = NULL) + 
  scale_x_continuous(expand = expansion(add = c(0.35, 0.5)), labels = NULL)

p3 <- plot_grid(p2, p1, align = "v",
                ncol = 1,
                rel_heights = c(2/10, 8/10),
                axis = "lr")

ggsave(
  "C:/Users/marce/Desktop/barplot_clut_dend.pdf",
  plot = p3,
  dpi = 10,
  device = "pdf",
  width = 15,
  height = 8,
)

########### Dendogram + barplot, no pattern

p4 <- ggplot(data = df_joined_filtered, aes(x = Sample, y=Abundance)) +
  geom_bar(aes(fill = Species),
                   position = "fill",
                   stat="identity",
                   show.legend = TRUE) +
  ggplot2::scale_fill_manual(values=colours_vec) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                 axis.text.y = ggplot2::element_text(size=12),
                 legend.text = element_text(size=12)) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")))

p4

p5 <- ggdendrogram(hc, rotate = 0,
                   leaf_labels = F) +
  theme(plot.margin = margin(t = 50,  # Top margin
                             r = 0,  # Right margin
                             b = -20,  # Bottom margin
                             l = 0)) + # Left margin
  scale_y_continuous(expand = expansion(add = c(1, 1)), labels = NULL) + 
  scale_x_continuous(expand = expansion(add = c(0.35, 0.5)), labels = NULL)

p6 <- plot_grid(p5, p4, align = "v",
                ncol = 1,
                rel_heights = c(2/10, 8/10),
                axis = "lr")

p6

ggsave(
  "C:/Users/marce/Desktop/barplot_clut_dend_no_pat.pdf",
  plot = p6,
  dpi = 10,
  device = "pdf",
  width = 15,
  height = 8,
)

