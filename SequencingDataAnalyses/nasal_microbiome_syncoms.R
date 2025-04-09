# Load scripts
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

### ... ###

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

# Convert OTU Table to strain level table.

strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA32", col_names = TRUE)

strain_ft <- merge_abundance_by_strain(otu_table_sctp_filt, strain_data)

# To use species-level data
otu_table <- otu_table_sctp_filt

# To use strain-level data
otu_table <- strain_ft

sc4 <- otu_table[1:12]
sc7 <- otu_table[13:24]
sc9 <- otu_table[25:36]
sc10 <- otu_table[37:48]
sc11 <- otu_table[49:60]
sc12 <- otu_table[61:72]
sc13 <- otu_table[73:84]
sc14 <- otu_table[85:96]
sc20 <- otu_table[97:108]
sc23 <- otu_table[109:120]
sc24 <- otu_table[121:132]
sc25 <- otu_table[133:144]
sc26 <- otu_table[145:156]
sc28 <- otu_table[157:168]
sc32 <- otu_table[169:180]
sc36 <- otu_table[181:192]
sc42 <- otu_table[193:204]
sc43 <- otu_table[205:216]
sc47 <- otu_table[217:228]
sc53 <- otu_table[229:240]

#barplot_w_strain_data2(sc4)
#barplot_w_strain_data2(sc7)
#barplot_w_strain_data2(sc9)

colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3", "springgreen4")

barplots_grid(feature_tables = list(sc4, sc7, sc9, sc10, sc11,
                                    sc12, sc13, sc14,sc20,sc23),
              strains = TRUE, shared_samples = FALSE,
              experiments_names = c("SC4", "SC7", "SC9", "SC10","SC11",
                                    "SC12", "SC13", "SC14", "SC20","SC23"),
              x_axis_title_size = 10, x_axis_text_size = 6,
              y_axis_title_size = 10, y_axis_text_size = 6,
              legend_pos = "bottom", legend_cols = 2,
              legend_title_size = 9, legend_text_size = 7,
              legend_key_size = 0.3, colour_palette = colours_vec)

#, "Strain 4" = "hexagonal", "Strain 5" = "pythagorean"
#, 0.2, 0.2
barplots_grid(feature_tables = list(sc24, sc25, sc26, sc28, sc32,
                                    sc36, sc42, sc43, sc47, sc53),
              strains = TRUE, shared_samples = FALSE,
              experiments_names = c("SC24", "SC25", "SC26", "SC28", "SC32",
                                    "SC36", "sc42", "SC43", "SC47","SC53"),
              x_axis_title_size = 10, x_axis_text_size = 6,
              y_axis_title_size = 10, y_axis_text_size = 6,
              legend_pos = "bottom", legend_cols = 2,
              legend_title_size = 9, legend_text_size = 7,
              legend_key_size = 0.3, colour_palette = colours_vec)





barplots_grid(feature_tables = list(sc4, sc7, sc9, sc10, sc11),
              strains = TRUE, shared_samples = FALSE,
              experiments_names = c("SC4", "SC7", "SC9", "SC10","SC11"),
              x_axis_title_size = 10, x_axis_text_size = 6,
              y_axis_title_size = 10, y_axis_text_size = 6,
              legend_pos = "bottom", legend_cols = 2,
              legend_title_size = 9, legend_text_size = 7,
              legend_key_size = 0.3, colour_palette = colours_vec)
