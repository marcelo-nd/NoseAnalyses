# Load helping functions
source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")

# Screening

otu_table_screening <- read.csv("F:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1)

colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8",
                                   "SC9", "SC10", "SC11", "SC12", "SC13", "SC14", "SC15",
                                   "SC16", "SC17", "SC18", "SC19", "SC20", "SC21", "SC22",
                                   "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29",
                                   "SC30", "SC31", "SC32", "SC33", "SC34", "SC35", "SC36",
                                   "SC37", "SC38", "SC39", "SC40", "SC41", "SC42", "SC43",
                                   "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50",
                                   "SC51", "SC52", "SC53")

barplot_from_feature_table(otu_table_screening)

# Remove SC39, contaminated with P. aueruginosa
df_filtered <- otu_table_screening[, -39]
# Remove SC34, contains Micrococcus luteus.
df_filtered <- df_filtered[, -34]
# Remove SC31. Apparent bad. S. aureus inoculation.
df_filtered <- df_filtered[, -31]
# Remove SC5, not correct species included
df_filtered <- df_filtered[, -5]
# Change name of C. keffirresidenti to C. cuberculostearicum
rownames(df_filtered)[rownames(df_filtered) == "Corynebacterium kefirresidentii"] <- "Corynebacterium tuberculostearicum"
# Remove Unassigned Readcounts
df_filtered <- df_filtered[-21,]

# Remove species with no counts

df_filtered <- filter_otus_by_counts_col_counts(df_filtered,min_count = 10, col_number = 1)

otu_table_adjusted <- df_filtered
otu_table_adjusted[7, 25:29] <- otu_table_adjusted[7, 25:29] / 5

barplot_from_feature_table(otu_table_adjusted)

###### Plots with fixed colour palette 
colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#290f76",
                "lightblue1","brown1", "olivedrab3", "darkorange3")

barplot_from_feature_table(feature_table = otu_table_adjusted, colour_palette = colours_vec)

barplot_from_feature_table_sorted(feature_table = otu_table_adjusted, colour_palette = colours_vec,
                                  sort_type = "species_abundance", species_to_sort = "Staphylococcus aureus")

barplot_from_feature_table_sorted(feature_table = otu_table_adjusted, colour_palette = colours_vec,
                                  sort_type = "species_abundance", species_to_sort = "Corynebacterium propinquum")

barplot_from_feature_table_sorted(feature_table = otu_table_adjusted, colour_palette = colours_vec,
                                  sort_type = "similarity")

############
##### Timepoints

otu_table_SC_tp <- read.csv("F:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv", row.names=1)

cf <- colnames(otu_table_SC_tp)

ordered_names <- cf[order(nchar(cf), cf)]

otu_table_SC_tp <- otu_table_SC_tp[, ordered_names]

colnames(otu_table_SC_tp) <- c("SC4_T1_R1", "SC4_T1_R2", "SC4_T1_R3",
                               "SC4_T2_R1", "SC4_T2_R2", "SC4_T2_R3",
                               "SC4_T3_R1", "SC4_T3_R2", "SC4_T3_R3",
                               "SC4_TF_R1", "SC4_TF_R2", "SC4_TF_R3",
                               "SC7_T1_R1", "SC7_T1_R2", "SC7_T1_R3",
                               "SC7_T2_R1", "SC7_T2_R2", "SC7_T2_R3",
                               "SC7_T3_R1", "SC7_T3_R2", "SC7_T3_R3",
                              "SC7_TF_R1", "SC7_TF_R2", "SC7_TF_R3",
                              "SC13_T1_R1", "SC13_T1_R2", "SC13_T1_R3",
                              "SC13_T2_R1", "SC13_T2_R2", "SC13_T2_R3",
                              "SC13_T3_R1", "SC13_T3_R2", "SC13_T3_R3",
                              "SC13_TF_R1", "SC13_TF_R2", "SC13_TF_R3",
                              "SC20_TF_R1", "SC20_TF_R2", "SC20_TF_R3",
                              "SC24_TF_R1", "SC24_TF_R2", "SC24_TF_R3",
                              "SC25_T1_R1", "SC25_T1_R2", "SC25_T1_R3",
                              "SC25_T2_R1", "SC25_T2_R2", "SC25_T2_R3",
                              "SC25_T3_R1", "SC25_T3_R2", "SC25_T3_R3",
                              "SC25_TF_R1", "SC25_TF_R2", "SC25_TF_R3",
                              "SC32_TF_R1", "SC32_TF_R2", "SC32_TF_R3",
                              "SC43_T1_R1", "SC43_T1_R2", "SC43_T1_R3",
                              "SC43_T2_R1", "SC43_T2_R2", "SC43_TF_R3",
                              "SC43_T3_R1", "SC43_T3_R2", "SC43_T2_R3",
                              "SC43_TF_R1", "SC43_TF_R2", "SC43_T3_R3",
                              "SC47_TF_R1", "SC47_TF_R2", "SC47_TF_R3",
                              "SC53_TF_R1", "SC53_TF_R2", "SC53_TF_R3",
                              "SC11_T1_R1", "SC11_T1_R2", "SC11_T1_R3",
                              "SC11_T2_R1", "SC11_T2_R2", "SC11_T2_R3",
                              "SC11_T3_R1", "SC11_T3_R2", "SC11_T3_R3",
                              "SC11_TF_R1", "SC11_TF_R2", "SC11_TF_R3",
                              "SC25_TF_R4"
                              )

otu_table_SC_tp_filt <- filter_otus_by_counts_col_counts(otu_table_SC_tp,
                                                             min_count = 10,
                                                             col_number = 1)

# Remove unnassigned reads
otu_table_SC_tp_filt <- otu_table_SC_tp_filt[-13,]

otu_table_SC_tp_filt["Corynebacterium tuberculostearicum", ] <- otu_table_SC_tp_filt["Corynebacterium tuberculostearicum", ] + otu_table_SC_tp_filt["Corynebacterium kefirresidentii", ]
otu_table_SC_tp_filt <- otu_table_SC_tp_filt[!rownames(otu_table_SC_tp_filt) %in% "Corynebacterium kefirresidentii", ]

otu_table_SC_tp_filt <- otu_table_SC_tp_filt[!rownames(otu_table_SC_tp_filt) %in% "Anaerococcus octavius", ]

# Change name of C. keffirresidenti to C. cuberculostearicum
#rownames(otu_table_SC_tp_filt)[rownames(otu_table_SC_tp_filt) == "Corynebacterium kefirresidentii"] <- "Corynebacterium tuberculostearicum"

barplot_from_feature_table(otu_table_SC_tp_filt)

colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                "lightblue1", "springgreen4", "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_table(otu_table_SC_tp_filt, colour_palette = colours_vec)

# Create subtables for each SC
sc4 <- otu_table_SC_tp_filt[1:12]
sc7 <- otu_table_SC_tp_filt[13:24]
sc13 <- otu_table_SC_tp_filt[25:36]
sc20 <- otu_table_SC_tp_filt[37:39]
sc24 <- otu_table_SC_tp_filt[40:42]
sc25 <- otu_table_SC_tp_filt[43:54]
sc32 <- otu_table_SC_tp_filt[55:57]
sc43 <- otu_table_SC_tp_filt[58:69]
sc47 <- otu_table_SC_tp_filt[70:72]
sc53 <- otu_table_SC_tp_filt[73:75]
sc11 <- otu_table_SC_tp_filt[76:87]

barplot_from_feature_tables(feature_tables = list(sc4, sc7, sc13,
                                                  sc20, sc24, sc25,
                                                  sc32, sc43, sc47, sc53, sc11),
                            experiments_names = c("SC4", "SC7", "SC13",
                                                  "SC20", "SC24", "SC25",
                                                  "SC32", "SC43", "SC47", "SC53", "sc11"),
                            colour_palette = colours_vec)

colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc4, sc25),
                            experiments_names = c("SC4", "SC25"),
                            colour_palette = colours_vec)

##### Time series barplots
# SC4
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

# SC7
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

# SC25
sc25_t1 <- sc25[1:3]
sc25_t2 <- sc25[4:6]
sc25_t3 <- sc25[7:9]
sc25_t4 <- sc25[10:12]

colours_vec <- c("#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc25_t1, sc25_t2, sc25_t3, sc25_t4),
                            experiments_names = c("SC25_T1", "SC25_T2", "SC25_T3", "SC25_T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC25",
                            x_axis_title_size = 10, y_axis_title_size = 10)

# SC43
sc43_t1 <- sc43[1:3]
sc43_t2 <- sc43[4:6]
sc43_t3 <- sc43[7:9]
sc43_t4 <- sc43[10:12]

colours_vec <- c("#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc43_t1, sc43_t2, sc43_t3, sc43_t4),
                            experiments_names = c("SC43_T1", "SC43_T2", "SC43_T3", "SC43_T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC43",
                            x_axis_title_size = 10, y_axis_title_size = 10)

# SC11
sc11_t1 <- sc11[1:3]
sc11_t2 <- sc11[4:6]
sc11_t3 <- sc11[7:9]
sc11_t4 <- sc11[10:12]

colours_vec <- c("#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3")

barplot_from_feature_tables(feature_tables = list(sc11_t1, sc11_t2, sc11_t3, sc11_t4),
                            experiments_names = c("SC11_T1", "SC11_T2", "SC11_T3", "SC11_T4"),
                            colour_palette = colours_vec, x_axis_text_size = 8, y_axis_text_size = 8,
                            legend_title_size = 10, legend_text_size = 8, plot_title = "SC11",
                            x_axis_title_size = 10, y_axis_title_size = 10)

