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
# Remove SC34, contains. Micrococcus luteus
df_filtered <- df_filtered[, -34]
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

colnames(otu_table_SC_tp) <- c("SC4_TF_R1", "SC4_TF_R2", "SC4_TF_R3",
                              "SC7_TF_R1", "SC7_TF_R2", "SC7_TF_R3",
                              "SC13_TF_R1", "SC13_TF_R2", "SC13_TF_R3",
                              "SC20_TF_R1", "SC20_TF_R2", "SC20_TF_R3",
                              "SC24_TF_R1", "SC24_TF_R2", "SC24_TF_R3",
                              "SC25_TF_R1", "SC25_TF_R2", "SC25_TF_R3",
                              "SC32_TF_R1", "SC32_TF_R2", "SC32_TF_R3",
                              "SC43_TF_R1", "SC43_TF_R2", "SC43_TF_R3",
                              "SC47_TF_R1", "SC47_TF_R2", "SC47_TF_R3",
                              "SC53_TF_R1", "SC53_TF_R2", "SC53_TF_R3")

otu_table_SC_tp_filt <- filter_otus_by_counts_col_counts(otu_table_SC_tp,
                                                             min_count = 10,
                                                             col_number = 1)

# Remove unnassigned reads
otu_table_SC_tp_filt <- otu_table_SC_tp_filt[-10,]

# Change name of C. keffirresidenti to C. cuberculostearicum
rownames(otu_table_SC_tp_filt)[rownames(otu_table_SC_tp_filt) == "Corynebacterium kefirresidentii"] <- "Corynebacterium tuberculostearicum"

barplot_from_feature_table(otu_table_SC_tp_filt)

colors_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                "lightblue1","brown1", "olivedrab3", "darkorange3")

barplot_from_feature_table(otu_table_SC_tp_filt, colour_palette = colors_vec)

# Create subtables for each SC
sc4 <- otu_table_SC_tp_filt[1:3]
sc7 <- otu_table_SC_tp_filt[4:6]
sc13 <- otu_table_SC_tp_filt[7:9]
sc20 <- otu_table_SC_tp_filt[10:12]
sc24 <- otu_table_SC_tp_filt[13:15]
sc25 <- otu_table_SC_tp_filt[16:18]
sc32 <- otu_table_SC_tp_filt[19:21]
sc43 <- otu_table_SC_tp_filt[22:24]
sc47 <- otu_table_SC_tp_filt[25:27]
sc53 <- otu_table_SC_tp_filt[28:30]

barplot_from_feature_tables(feature_tables = list(sc4, sc7, sc13,
                                                  sc20, sc24, sc25,
                                                  sc32, sc43, sc47, sc53),
                            experiments_names = c("SC4", "SC7", "SC13",
                                                  "SC20", "SC24", "SC25",
                                                  "SC32", "SC43", "SC47", "SC53"),
                            colour_palette = colours_vec)

