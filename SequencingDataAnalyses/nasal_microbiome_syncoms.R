# Load scripts
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_graphing.R")

##### Screening Results
otu_table_screening <- read.csv("E:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1)


colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC17", "SC18", "SC19", "SC20", "SC21",
                                   "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30", "SC31",
                                   "SC32", "SC33", "SC35", "SC36", "SC37", "SC38", "SC40", "SC41", "SC42", "SC43",
                                   "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50", "SC51", "SC52", "SC53")

# List of species to remove (they did not grow in any of the SynComs, and Unassigned reads)
species_to_remove <- c("Anaerococcus octavius", "Cutibacterium acnes", "Unassigned")

ot_scree_filtered <- remove_feature_by_prefix(otu_table_screening, species_to_remove)

write.csv(x = ot_scree_filtered,
          file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/scree_ot_filtered.csv",
          row.names = T,
          quote = F)

colours_vec <- c("#ffe599", "dodgerblue4", "blueviolet", "#CC79A7","mediumspringgreen",
                 "lightblue1","#EF5B5B", "olivedrab3", "#e89d56")

##### Screening Barplot
screening_barplot <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "none", colour_palette = colours_vec, legend_pos = "bottom", legend_cols = 3)
screening_barplot

##### Screening Barplot sorted by S. aureus abundance
screening_barplot_sau_sorted <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "feature_value", feature_to_sort = "Staphylococcus aureus",
                                                colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)
screening_barplot_sau_sorted

##### Screening Barplot sorted by Similarity
screening_barplot_sim_sorted <- barplot_from_feature_table(feature_table = ot_scree_filtered, sort_type = "similarity",
                                                colour_palette = colours_vec, legend_pos = "bottom", legend_cols = 4)
screening_barplot_sim_sorted

ggsave(
  "C:/Users/marce/Desktop/bar_plot_scree.pdf",
  plot = screening_barplot,
  dpi = 10, device = "pdf", width = 15, height = 6
)

##### Add dendograms to barplots
dendo_plot <- dendrogram_from_feature_table(ot_scree_filtered)

barplor_w_dendo <- cowplot::plot_grid(dendo_plot, screening_barplot_sim_sorted, align = "v",
                         ncol = 1,
                         rel_heights = c(2/10, 8/10),
                         axis = "lr")

barplor_w_dendo


# Screening barplot with strain information
# Convert OTU Table to strain level table.

strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA32", col_names = TRUE)

strain_ft <- merge_abundance_by_strain(ot_scree_filtered, strain_data)

screening_barplot <- barplot_from_feature_table(feature_table = strain_ft, sort_type = "none", strains = TRUE,
                                                colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)

screening_barplot


