library(dplyr)
library(ggplot2)
library(readxl)

############## Real run
source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")

# Read table, todo order or not
fia_table <- read_fia_table("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/20250115_FIA-Data_lc_SCs_exp.xlsx")
fia_table <- fia_table[order(row.names(fia_table)), ]
write.csv(x = fia_table, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/feature_table.csv")

LC_metadata <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/Metabolomics tubes annotation - Nov 2024_2.xlsx", sheet = "Tabelle1")
LC_metadata <- LC_metadata[order(LC_metadata$sample), ]
write.csv(x = LC_metadata, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/metadata.csv")

# Extract Interesting features
ft_subset <- extract_features_comparison(feature_table = fia_table, sig_features = c("L-Alanine[M+H]+", "L-Arginine[M+H]+", "L-Asparagine[M+H]+",
                                                                                     "L-Aspartic acid[M+H]+", "L-Cystine[M+H]+", "L-Glutamic acid[M+H]+",
                                                                                     "L-Glutamine[M+H]+", "L-Histidine[M+H]+", "L-Isoleucine[M+H]+",
                                                                                     "L-Leucine[M+H]+", "L-Lysine[M+H]+", "L-Methionine[M+H]+",
                                                                                     "L-Phenylalanine[M+H]+", "L-Proline[M+H]+", "L-Serine[M+H]+",
                                                                                     "L-Threonine[M+H]+", "L-Tryptophan[M+H]+", "L-Tyrosine[M+H]+",
                                                                                     "L-Valine[M+H]+"), columns_to_preserve = NULL)
write.csv(x = ft_subset, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/ft_subset.csv")

# Normalize
ft_od <- normalize_by_od(feature_table = ft_subset, metadata_table = LC_metadata, select_by = "all", samples_group_to_exclude = "PBS")

# Log transform
ft_log <- log2_convert(ft_subset)

# Normalize to treatment
#ft_log_norm <- normalize_table_to_treatment(metabolites_table = ft_log, metadata_table = LC_metadata, samples_norm = colnames(select(ft_log, contains("PBS"))))
ft_log_norm <- normalize_table_to_treatment(metabolites_table = ft_log, metadata_table = LC_metadata, samples_norm = "PBS")

write.csv(x = ft_log_norm, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/ft_log_norm.csv")
# Graph
graph_metabolites(ft_log_norm, y1 = -8, y2 = 8, dotsize = 25, binwidth = 0.01)

# DO PCA
ft_pca(fia_table, metadata_table = LC_metadata, grouping_col = syncom)



































































feature_table <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_Test_230824/20240822_FIA-Data_LC_Test_pos.xlsx")

#feature_table <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_Test_230824/20240822_MergePos_LC_Test.xlsx")

LC_metadata <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_Test_230824/LC_metadata.xlsx")

transdata <- t(feature_table)
transdata <- as.data.frame(transdata)
colnames(transdata) <- transdata[2,]
transdata <- transdata[5:nrow(transdata), ]

#fiadata <- select(X20240822_FIA_Data_LC_Test_pos, c(2, 5:ncol(X20240822_FIA_Data_LC_Test_pos)))
#metNames <- fiadata$MetName
#fiadata <- fiadata[2:ncol(fiadata)]
#row.names(fiadata) <- metNames

new_sample_names <- paste(LC_metadata$syncom, LC_metadata$group, LC_metadata$replicate, sep = "_")

in_feat_table <- extract_features_comparison(feature_table = transdata, sig_features = c("L-Alanine[M+H]+", "L-Arginine[M+H]+", "L-Asparagine[M+H]+",
                                                                                          "L-Aspartic acid[M+H]+", "L-Cystine[M+H]+", "L-Glutamic acid[M+H]+",
                                                                                         "L-Glutamine[M+H]+", "L-Histidine[M+H]+", "L-Isoleucine[M+H]+",
                                                                                         "L-Leucine[M+H]+", "L-Lysine[M+H]+", "L-Methionine[M+H]+",
                                                                                         "L-Phenylalanine[M+H]+", "L-Proline[M+H]+", "L-Serine[M+H]+",
                                                                                         "L-Threonine[M+H]+", "L-Tryptophan[M+H]+", "L-Tyrosine[M+H]+",
                                                                                         "L-Valine[M+H]+"), columns_to_preserve = NULL)
in_feat_table <- in_feat_table

colnames(in_feat_table) <- c(new_sample_names, "Metabolite")


metab_flooding <- select(in_feat_table, contains("flooding"), "Metabolite")

metab_flooding_log <- log2_convert(metab_flooding)

norm_table_metab_flooding <- normalize_table_to_treatment(metabolites_table = metab_flooding_log, samples_norm = colnames(select(metab_flooding, contains("Control"))))

metab_flooding_g <- tidyr::gather(data = norm_table_metab_flooding, key = "Sample", value = "Change", colnames(norm_table_metab_flooding[,2:ncol(norm_table_metab_flooding)]))

metab_flooding_g$sample_type <- sapply(strsplit(metab_flooding_g$Sample, "_"), "[", 1)

graph_metabolites(metab_flooding_g, , y1 = -7.5, y2 = 3, dotsize = 25, binwidth = 0.01)



metab_scrapping <- select(in_feat_table, contains("scrapping"), "Metabolite")

metab_scrapping_log <- log2_convert(metab_scrapping)

norm_table_metab_scrapping <- normalize_table_to_treatment(metabolites_table = metab_scrapping_log, samples_norm = colnames(select(metab_scrapping, contains("Control"))))

metab_scrapping_g <- tidyr::gather(data = norm_table_metab_scrapping, key = "Sample", value = "Change", colnames(norm_table_metab_scrapping[,2:ncol(norm_table_metab_scrapping)]))

metab_scrapping_g$sample_type <- sapply(strsplit(metab_scrapping_g$Sample, "_"), "[", 1)

graph_metabolites(metab_scrapping_g, , y1 = -3, y2 = 7.5, dotsize = 25, binwidth = 0.01)


























# Data forming
ft_log_norm_g <- tidyr::gather(data = ft_log_norm, key = "Sample", value = "Change", colnames(ft_log_norm[,2:ncol(ft_log_norm)]))

ft_log_norm_g$sample_type <- sapply(strsplit(ft_log_norm_g$Sample, ".", fixed = TRUE), "[", 1)

# Graphing
graph_metabolites(ft_log_norm_g, , y1 = -7.5, y2 = 3, dotsize = 25, binwidth = 0.01)





feature_table_lc <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/20250115_FIA-Data_lc_SCs_exp.xlsx")



LC_metadata <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/Metabolomics tubes annotation - Nov 2024_2.xlsx", sheet = "Tabelle1")

LC_metadata <- LC_metadata[order(LC_metadata$sample), ]

transdata <- t(feature_table_lc)
colnames(transdata) <- transdata[2,]
transdata <- transdata[5:nrow(transdata), ]
transdata <- transdata[order(row.names(transdata)), ]
transdata <- as.data.frame(transdata)

#transdata2 <- cbind(LC_metadata["od"], transdata)

#df_normalized <- transdata2 %>% mutate(across(where(is.numeric), ~ . / od))

# Normalize by OD
transdata <- extract_features_comparison(feature_table = transdata, sig_features = c("L-Alanine[M+H]+", "L-Arginine[M+H]+", "L-Asparagine[M+H]+",
                                                                                         "L-Aspartic acid[M+H]+", "L-Cystine[M+H]+", "L-Glutamic acid[M+H]+",
                                                                                         "L-Glutamine[M+H]+", "L-Histidine[M+H]+", "L-Isoleucine[M+H]+",
                                                                                         "L-Leucine[M+H]+", "L-Lysine[M+H]+", "L-Methionine[M+H]+",
                                                                                         "L-Phenylalanine[M+H]+", "L-Proline[M+H]+", "L-Serine[M+H]+",
                                                                                         "L-Threonine[M+H]+", "L-Tryptophan[M+H]+", "L-Tyrosine[M+H]+",
                                                                                         "L-Valine[M+H]+"), columns_to_preserve = NULL)

transdata2 <- t(transdata[1:ncol(transdata) - 1])

transdata2 <- as.data.frame(transdata2)

ft_od <- normalize_by_od(feature_table = transdata2, metadata_table = LC_metadata, od_column = "od", select_by = "all")

ft_log <- log2_convert(ft_od)

ft_log_norm <- normalize_table_to_treatment(metabolites_table = ft_log, samples_norm = colnames(select(ft_log, contains("PBS"))))

ft_log_norm_g <- tidyr::gather(data = ft_log_norm, key = "Sample", value = "Change", colnames(ft_log_norm[,2:ncol(ft_log_norm)]))

ft_log_norm_g$sample_type <- sapply(strsplit(ft_log_norm_g$Sample, ".", fixed = TRUE), "[", 1)

graph_metabolites(ft_log_norm_g, , y1 = -7.5, y2 = 3, dotsize = 25, binwidth = 0.01)
