# Load helping functions
source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")

###################################################################################################
########################################## Screening ##############################################
otu_table_screening <- read.csv("E:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1)

colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10",
                                   "SC11", "SC12", "SC13", "SC14", "SC15", "SC17", "SC18", "SC19", "SC20", "SC21",
                                   "SC22", "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29", "SC30", "SC31",
                                   "SC32", "SC33", "SC35", "SC36", "SC37", "SC38", "SC40", "SC41", "SC42", "SC43",
                                   "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50", "SC51", "SC52", "SC53")

# Remove Unassigned Readcounts
ot_scree_filtered <- otu_table_screening[-10,]

# Remove species with no counts
ot_scree_filtered <- filter_otus_by_counts_col_counts(ot_scree_filtered,min_count = 10, col_number = 1)

# Set colour palette 
colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#290f76",
                "lightblue1","brown1", "olivedrab3", "darkorange3")

##### Screening Barplot
barplot_from_feature_table(feature_table = ot_scree_filtered, colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)

ggsave(
  "C:/Users/marce/Desktop/bar_plot_scree.pdf",
  plot = barplot_from_feature_table(feature_table = ot_scree_filtered, colour_palette = colours_vec, legend_pos = "right", legend_cols = 1),
  dpi = 10,
  device = "pdf",
  width = 15,
  height = 6
)

##### Screening Barplots sorted by S. aureus abundance and Similarity

barplot_from_feature_table_sorted(feature_table = ot_scree_filtered, sort_type = "species_abundance", species_to_sort = "Staphylococcus aureus",
                                  colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)

barplot_from_feature_table_sorted(feature_table = ot_scree_filtered, sort_type = "similarity",
                                  colour_palette = colours_vec, legend_pos = "right", legend_cols = 1)

write.csv(x = ot_scree_filtered,
          file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/scree_ot_filtered.csv",
          row.names = T,
          quote = F)

##### Screening Barplots with strain information.


####################################################################################################
########################################## Timepoints ##############################################
# Read table with timepoint data
otu_table_sctp <- read.csv("E:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv", row.names=1)

otu_table_sctp_sorted <- sort_nanopore_table_by_barcodes(df = otu_table_sctp,
                                                       new_names = c("SC4_T1_R1", "SC4_T1_R2", "SC4_T1_R3",
                                                                     "SC4_T2_R1", "SC4_T2_R2", "SC4_T2_R3",
                                                                     "SC4_T3_R1", "SC4_T3_R2", "SC4_T3_R3",
                                                                     "SC4_TF_R1", "SC4_TF_R2", "SC4_TF_R3",
                                                                     "SC7_T1_R1", "SC7_T1_R2", "SC7_T1_R3",
                                                                     "SC7_T2_R1", "SC7_T2_R2", "SC7_T2_R3",
                                                                     "SC7_T3_R1", "SC7_T3_R2", "SC7_T3_R3",
                                                                     "SC7_TF_R1", "SC7_TF_R2", "SC7_TF_R3",
                                                                     "SC9_T1_R1", "SC9_T1_R2", "SC9_T1_R3",
                                                                     "SC9_T2_R1", "SC9_T2_R2", "SC9_T2_R3",
                                                                     "SC9_T3_R1", "SC9_T3_R2", "SC9_T3_R3",
                                                                     "SC9_TF_R1", "SC9_TF_R2", "SC9_TF_R3",
                                                                     "SC10_T1_R1", "SC10_T1_R2", "SC10_T1_R3",
                                                                     "SC10_T2_R1", "SC10_T2_R2", "SC10_T2_R3",
                                                                     "SC10_T3_R1", "SC10_T3_R2", "SC10_T3_R3",
                                                                     "SC10_TF_R1", "SC10_TF_R2", "SC10_TF_R3",
                                                                     "SC11_T1_R1", "SC11_T1_R2", "SC11_T1_R3",
                                                                     "SC11_T2_R1", "SC11_T2_R2", "SC11_T2_R3",
                                                                     "SC11_T3_R1", "SC11_T3_R2", "SC11_T3_R3",
                                                                     "SC11_TF_R1", "SC11_TF_R2", "SC11_TF_R3",
                                                                     "SC12_T1_R1", "SC12_T1_R2", "SC12_T1_R3",
                                                                     "SC12_T2_R1", "SC12_T2_R2", "SC12_T2_R3",
                                                                     "SC12_T3_R1", "SC12_T3_R2", "SC12_T3_R3",
                                                                     "SC12_TF_R1", "SC12_TF_R2", "SC12_TF_R3",
                                                                     "SC13_T1_R1", "SC13_T1_R2", "SC13_T1_R3",
                                                                     "SC13_T2_R1", "SC13_T2_R2", "SC13_T2_R3",
                                                                     "SC13_T3_R1", "SC13_T3_R2", "SC13_T3_R3",
                                                                     "SC13_TF_R1", "SC13_TF_R2", "SC13_TF_R3",
                                                                     "SC14_T1_R1", "SC14_T1_R2", "SC14_T1_R3",
                                                                     "SC14_T2_R1", "SC14_T2_R2", "SC14_T2_R3",
                                                                     "SC14_T3_R1", "SC14_T3_R2", "SC14_T3_R3",
                                                                     "SC14_TF_R1", "SC14_TF_R2", "SC14_TF_R3",
                                                                     "SC20_T1_R1", "SC20_T1_R2", "SC20_T1_R3",
                                                                     "SC20_T2_R1", "SC20_T2_R2", "SC20_T2_R3",
                                                                     "SC20_T3_R1", "SC20_T3_R2", "SC20_T3_R3",
                                                                     "SC20_TF_R1", "SC20_TF_R2", "SC20_TF_R3",
                                                                     "SC23_T1_R1", "SC23_T1_R2", "SC23_T1_R3",
                                                                     "SC23_T2_R1", "SC23_T2_R2", "SC23_T2_R3",
                                                                     "SC23_T3_R1", "SC23_T3_R2", "SC23_T3_R3",
                                                                     "SC23_TF_R1", "SC23_TF_R2", "SC23_TF_R3",
                                                                     "SC24_T1_R1", "SC24_T1_R2", "SC24_T1_R3",
                                                                     "SC24_T2_R1", "SC24_T2_R2", "SC24_T2_R3",
                                                                     "SC24_T3_R1", "SC24_T3_R2", "SC24_T3_R3",
                                                                     "SC24_TF_R1", "SC24_TF_R2", "SC24_TF_R3",
                                                                     "SC25_T1_R1", "SC25_T1_R2", "SC25_T1_R3",
                                                                     "SC25_T2_R1", "SC25_T2_R2", "SC25_T2_R3",
                                                                     "SC25_T3_R1", "SC25_T3_R2", "SC25_T3_R3",
                                                                     "SC25_TF_R1", "SC25_TF_R2", "SC25_TF_R3",
                                                                     "SC26_T1_R1", "SC26_T1_R2", "SC26_T1_R3",
                                                                     "SC26_T2_R1", "SC26_T2_R2", "SC26_T2_R3",
                                                                     "SC26_T3_R1", "SC26_T3_R2", "SC26_T3_R3",
                                                                     "SC26_TF_R1", "SC26_TF_R2", "SC26_TF_R3",
                                                                     "SC28_T1_R1", "SC28_T1_R2", "SC28_T1_R3",
                                                                     "SC28_T2_R1", "SC28_T2_R2", "SC28_T2_R3",
                                                                     "SC28_T3_R1", "SC28_T3_R2", "SC28_T3_R3",
                                                                     "SC28_TF_R1", "SC28_TF_R2", "SC28_TF_R3",
                                                                     "SC32_T1_R1", "SC32_T1_R2", "SC32_T1_R3",
                                                                     "SC32_T2_R1", "SC32_T2_R2", "SC32_T2_R3",
                                                                     "SC32_T3_R1", "SC32_T3_R2", "SC32_T3_R3",
                                                                     "SC32_TF_R1", "SC32_TF_R2", "SC32_TF_R3",
                                                                     "SC36_T1_R1", "SC36_T1_R2", "SC36_T1_R3",
                                                                     "SC36_T2_R1", "SC36_T2_R2", "SC36_T2_R3",
                                                                     "SC36_T3_R1", "SC36_T3_R2", "SC36_T3_R3",
                                                                     "SC36_TF_R1", "SC36_TF_R2", "SC36_TF_R3",
                                                                     "SC42_T1_R1", "SC42_T1_R2", "SC42_T1_R3",
                                                                     "SC42_T2_R1", "SC42_T2_R2", "SC42_T2_R3",
                                                                     "SC42_T3_R1", "SC42_T3_R2", "SC42_T3_R3",
                                                                     "SC42_TF_R1", "SC42_TF_R2", "SC42_TF_R3",
                                                                     "SC43_T1_R1", "SC43_T1_R2", "SC43_T1_R3",
                                                                     "SC43_T2_R1", "SC43_T2_R2", "SC43_T2_R3",
                                                                     "SC43_T3_R1", "SC43_T3_R2", "SC43_T3_R3",
                                                                     "SC43_TF_R1", "SC43_TF_R2", "SC43_TF_R3",
                                                                     "SC47_T1_R1", "SC47_T1_R2", "SC47_T1_R3",
                                                                     "SC47_T2_R1", "SC47_T2_R2", "SC47_T2_R3",
                                                                     "SC47_T3_R1", "SC47_T3_R2", "SC47_T3_R3",
                                                                     "SC47_TF_R1", "SC47_TF_R2", "SC47_TF_R3",
                                                                     "SC53_T1_R1", "SC53_T1_R2", "SC53_T1_R3",
                                                                     "SC53_T2_R1", "SC53_T2_R2", "SC53_T2_R3",
                                                                     "SC53_T3_R1", "SC53_T3_R2", "SC53_T3_R3",
                                                                     "SC53_TF_R1", "SC53_TF_R2", "SC53_TF_R3"))

# Remove species with no counts
otu_table_sctp_filt <- filter_otus_by_counts_col_counts(otu_table_sctp_sorted,
                                                             min_count = 10,
                                                             col_number = 1)

# Remove Anaerococcus octavius, it did not grow on any of the SCs
otu_table_sctp_filt <- otu_table_sctp_filt[!rownames(otu_table_sctp_filt) %in% "Anaerococcus octavius", ]

# Remove Cutibacterium acnes, it did not grow on any of the SCs
otu_table_sctp_filt <- otu_table_sctp_filt[!rownames(otu_table_sctp_filt) %in% "Cutibacterium acnes", ]

# Remove unnassigned reads
otu_table_sctp_filt <- otu_table_sctp_filt[-10,]

# Barplot with fixed palette
colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3", "springgreen4")

barplot_from_feature_table(otu_table_sctp_filt, colour_palette = colours_vec, legend_pos = "bottom",
                           legend_cols = 3, x_axis_text_size = 4, y_axis_title_size = 8, y_axis_text_size = 5)

write.table(x = otu_table_sctp_filt,
          file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/sc100_tp_ot_filtered.csv",
          row.names = TRUE, col.names = TRUE, quote = F, sep = ";")

##### Get relative abundance table. Useful for other analyses.
otu_table_sctp_filt_rel_ab <- calculate_relative_abundance(otu_table_sctp_filt)

write.table(x = otu_table_sctp_filt_rel_ab,
          file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/champs_ot_filtered_rel_ab.csv",
          row.names = T,
          quote = F)

# Create subtables for each SC
sc4 <- otu_table_sctp_filt[1:12]
sc7 <- otu_table_sctp_filt[13:24]
sc9 <- otu_table_sctp_filt[25:36]
sc10 <- otu_table_sctp_filt[37:48]
sc11 <- otu_table_sctp_filt[49:60]
sc12 <- otu_table_sctp_filt[61:72]
sc13 <- otu_table_sctp_filt[73:84]
sc14 <- otu_table_sctp_filt[85:96]
sc20 <- otu_table_sctp_filt[97:108]
sc23 <- otu_table_sctp_filt[109:120]
sc24 <- otu_table_sctp_filt[121:132]
sc25 <- otu_table_sctp_filt[133:144]
sc26 <- otu_table_sctp_filt[145:156]
sc28 <- otu_table_sctp_filt[157:168]
sc32 <- otu_table_sctp_filt[169:180]
sc36 <- otu_table_sctp_filt[181:192]
sc42 <- otu_table_sctp_filt[193:204]
sc43 <- otu_table_sctp_filt[205:216]
sc47 <- otu_table_sctp_filt[217:228]
sc53 <- otu_table_sctp_filt[229:240]

barplot_from_feature_tables(feature_tables = list(sc4, sc7, sc9, sc10, sc11,
                                                  sc12, sc13, sc14,sc20,sc23),
                            experiments_names = c("SC4", "SC7", "SC9", "SC10","sc11",
                                                  "SC12", "SC13", "SC14", "SC20","SC23"),
                            colour_palette = colours_vec,
                            x_axis_title_size = 9, x_axis_text_size = 5,
                            y_axis_title_size = 9, y_axis_text_size = 5,
                            legend_pos = "bottom", legend_cols = 2,
                            legend_title_size = 9, legend_text_size = 7,
                            legend_key_size = 0.3)

barplot_from_feature_tables(feature_tables = list(sc24, sc25, sc26, sc28, sc32,
                                                  sc36, sc42, sc43, sc47, sc53),
                            experiments_names = c("SC24", "SC25", "SC26", "SC28", "SC32",
                                                  "SC36", "sc42", "SC43", "SC47","SC53"),
                            colour_palette = colours_vec,
                            x_axis_title_size = 9, x_axis_text_size = 5,
                            y_axis_title_size = 9, y_axis_text_size = 5,
                            legend_pos = "bottom", legend_cols = 2,
                            legend_title_size = 9, legend_text_size = 7,
                            legend_key_size = 0.3)

barplot_with_replicates(sc4)
barplot_with_replicates(sc7)
barplot_with_replicates(sc9)
barplot_with_replicates(sc10)
barplot_with_replicates(sc11)
barplot_with_replicates(sc12)
barplot_with_replicates(sc13)
barplot_with_replicates(sc14)
barplot_with_replicates(sc20)
barplot_with_replicates(sc23)
barplot_with_replicates(sc24)
barplot_with_replicates(sc25)
barplot_with_replicates(sc26)
barplot_with_replicates(sc28)
barplot_with_replicates(sc32)
barplot_with_replicates(sc36)
barplot_with_replicates(sc42)
barplot_with_replicates(sc43)
barplot_with_replicates(sc47)
barplot_with_replicates(sc53)

##################################################################################################################################
################################################# Strain Level Barplots ##########################################################
library(ggpattern)

library(dendextend)

library(ggdendro)

library(cowplot)

barplot_w_strain_data()


##################################################################################################################################
# Read metadata
sc100_tp_metadata <- read.table(file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/sc100_tp_ot_filtered_metadata.csv",
                                header = TRUE,sep = ";", row.names = 1)


ft_pca_2(t(otu_table_sctp_filt_rel_ab), metadata_table = sc100_tp_metadata, grouping_col = "syncom")


feature_table <- read_ft_1("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/annotated_QuantTable_scaled.csv",
                           sort_by_names = TRUE)

feature_table2 <- read_ft_1("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_scaled.csv",
                           sort_by_names = TRUE)

rownames(feature_table2) <- gsub("\\.mzML$", "", rownames(feature_table2))


# Heatmaps
feature_table_heatmap(otu_table_sctp_filt_rel_ab, otu_table_sctp_filt_rel_ab)

feature_table_heatmap(otu_table_sctp_filt_rel_ab, t(feature_table2))

feature_table_heatmap(abundance_per_strain, abundance_per_strain)

feature_table_heatmap(abundance_per_strain, t(feature_table2))

# TO DO
feature_table_heatmap_w_sig(otu_table_sctp_filt_rel_ab, otu_table_sctp_filt_rel_ab)

feature_table_heatmap_w_sig(otu_table_sctp_filt_rel_ab, t(feature_table2))

# Random forest for pairs?
########################################################



##### Time series barplots
ot_scree_filtered <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/SynCom100/Results/scree_ot_filtered.csv",
                              header = T,
                              row.names=1)

nose_strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA32", col_names = TRUE)

nose_strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA31", col_names = TRUE)

test1 <- barplot_w_strain_data(otu_table = ot_scree_filtered, strain_data = nose_strain_data)

#barplot_w_strain_data(otu_table = otu_table_sctp_filt)
barplot_w_strain_data(otu_table = ot_scree_filtered, strain_data = nose_strain_data)

barplot_w_strain_data(otu_table = sc4_t1)

# Theoretical values
# Collapse the information by species
df_collapsed_3 <- strain_data %>%
  mutate(Species = sapply(strsplit(Species, " "), function(x) paste(x[1:2], collapse = " "))) %>% # Extract species name
  group_by(Species) %>%
  summarise(across(starts_with("SC"), max)) %>% # Take max per sample to represent strain
  ungroup()

spps <- df_collapsed_3$Species

df_collapsed_3 <- select(df_collapsed_3, -1)

rownames(df_collapsed_3) <- spps

# SC4
sc4_t0 <- select(df_collapsed_3, "SC4")

sc4_t0 <- as.data.frame(sc4_t0)

rownames(sc4_t0) <- spps

sc4_t0 <- filter_species_by_col_counts(sc4_t0, min_count = 1, col_number = 1)

sc4_t1 <- sc4[1:3]
sc4_t2 <- sc4[4:6]
sc4_t3 <- sc4[7:9]
sc4_t4 <- sc4[10:12]

colours_vec <- c("gold3", "#053f73", "blueviolet","#6279B8",
                 "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc4_t1, sc4_t2, sc4_t3, sc4_t4),
                            experiments_names = c("SC4_T1", "SC4_T2", "SC4_T3", "SC4_T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC4",
                            x_axis_title_size = 10, y_axis_title_size = 10)

colours_vec <- c("#4B1E19", "gold3", "#053f73", "blueviolet", "#CC79A7", "#14471E", "#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3", "springgreen4")

barplot_from_feature_tables(feature_tables = list(sc4_t0, sc4_t1, sc4_t2, sc4_t3, sc4_t4),
                            experiments_names = c("T0", "T1", "T2", "T3", "T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC4",
                            x_axis_title_size = 10, y_axis_title_size = 10, legend_pos = "bottom", legend_cols = 4)

# SC7
sc7_t0 <- select(df_collapsed_3, "SC7")

sc7_t0 <- as.data.frame(sc7_t0)

rownames(sc7_t0) <- spps

sc7_t0 <- filter_species_by_col_counts(sc7_t0, min_count = 1, col_number = 1)

sc7_t1 <- sc7[1:3]
sc7_t2 <- sc7[4:6]
sc7_t3 <- sc7[7:9]
sc7_t4 <- sc7[10:12]

colours_vec <- c("#053f73", "blueviolet","#6279B8",
                 "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc7_t1, sc7_t2, sc7_t3, sc7_t4),
                            experiments_names = c("SC7_T1", "SC7_T2", "SC7_T3", "SC7_T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC7",
                            x_axis_title_size = 10, y_axis_title_size = 10)

colours_vec <- c("#4B1E19", "gold3", "#053f73", "blueviolet", "#CC79A7", "#14471E", "#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3", "springgreen4")

barplot_from_feature_tables(feature_tables = list(sc7_t0, sc7_t1, sc7_t2, sc7_t3, sc7_t4),
                            experiments_names = c("T0", "T1", "T2", "T3", "T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC7",
                            x_axis_title_size = 10, y_axis_title_size = 10, legend_pos = "bottom", legend_cols = 4)

# SC11
sc11_t0 <- select(df_collapsed_3, "SC11")

sc11_t0 <- as.data.frame(sc11_t0)

rownames(sc11_t0) <- spps

sc11_t0 <- filter_species_by_col_counts(sc11_t0, min_count = 1, col_number = 1)


sc11_t1 <- sc11[1:3]
sc11_t2 <- sc11[4:6]
sc11_t3 <- sc11[7:9]
sc11_t4 <- sc11[10:12]

colours_vec <- c("gold3", "#053f73", "blueviolet", "#6279B8",
                 "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc11_t1, sc11_t2, sc11_t3, sc11_t4),
                            experiments_names = c("SC11_T1", "SC11_T2", "SC11_T3", "SC11_T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC11",
                            x_axis_title_size = 10, y_axis_title_size = 10)

colours_vec <- c("#4B1E19", "gold3", "#053f73", "blueviolet", "#CC79A7", "#14471E", "#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3", "springgreen4", "black")

barplot_from_feature_tables(feature_tables = list(sc11_t0, sc11_t1, sc11_t2, sc11_t3, sc11_t4),
                            experiments_names = c("T0", "T1", "T2", "T3", "T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC11",
                            x_axis_title_size = 10, y_axis_title_size = 10, legend_pos = "bottom", legend_cols = 5)

# SC25
sc25_t0 <- select(df_collapsed_3, "SC25")

sc25_t0 <- as.data.frame(sc25_t0)

rownames(sc25_t0) <- spps

sc25_t0 <- filter_species_by_col_counts(sc25_t0, min_count = 1, col_number = 1)

sc25_t1 <- sc25[1:3]
sc25_t2 <- sc25[4:6]
sc25_t3 <- sc25[7:9]
sc25_t4 <- sc25[10:12]

colours_vec <- c("#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc25_t1, sc25_t2, sc25_t3, sc25_t4),
                            experiments_names = c("T1", "T2", "SC25_T3", "SC25_T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC25",
                            x_axis_title_size = 10, y_axis_title_size = 10)

colours_vec <- c("#053f73", "blueviolet", "#CC79A7", "#14471E", "#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3", "springgreen4")

barplot_from_feature_tables(feature_tables = list(sc25_t0, sc25_t1, sc25_t2, sc25_t3, sc25_t4),
                            experiments_names = c("T0", "T1", "T2", "T3", "T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC25",
                            x_axis_title_size = 10, y_axis_title_size = 10, legend_pos = "bottom", legend_cols = 4)

# SC43
sc43_t0 <- select(df_collapsed_3, "SC43")
sc43_t0 <- as.data.frame(sc43_t0)
rownames(sc43_t0) <- spps
sc43_t0 <- sc43_t0[-3, , drop = F]

sc43_t0 <- filter_species_by_col_counts(sc43_t0, min_count = 1, col_number = 1)

sc43_t1 <- sc43[1:3]
sc43_t2 <- sc43[4:6]
sc43_t3 <- sc43[7:9]
sc43_t4 <- sc43[10:12]
sc43_t4 <- sc43_t4[-2,]

colours_vec <- c("gold3", "blueviolet", "#6279B8",
                 "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc43_t1, sc43_t2, sc43_t3, sc43_t4),
                            experiments_names = c("T1", "T2", "T3", "T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC43",
                            x_axis_title_size = 10, y_axis_title_size = 10, legend_cols = 2)

colours_vec <- c("#4B1E19", "gold3", "#14471E", "#6279B8",
                  "brown1", "olivedrab3", "darkorange3", "springgreen4")

barplot_from_feature_tables(feature_tables = list(sc43_t0, sc43_t1, sc43_t2, sc43_t3, sc43_t4),
                            experiments_names = c("T0", "T1", "T2", "T3", "T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC43",
                            x_axis_title_size = 10, y_axis_title_size = 10, legend_pos = "bottom", legend_cols = 3)
  