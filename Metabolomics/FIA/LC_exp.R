source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")

# Read table
fia_table <- read_fia_table("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/20250115_FIA-Data_lc_SCs_exp.xlsx",
                            sort_table = TRUE, fix_names = FALSE)

write.csv(x = fia_table, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/feature_table.csv")

# Read metadata
LC_metadata <- read_metadata_xls(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/Metabolomics tubes annotation - Nov 2024_2.xlsx",
                                 sheet = "Tabelle1", sort_table = TRUE)

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

##### OD normalized + log transformed + treatment-based log fold change
# Normalize by OD
ft_od <- normalize_by_od(feature_table = ft_subset, metadata_table = LC_metadata, grouping_variable = "syncom", samples_group_to_exclude = "PBS")

# Log transform 
ft_od_log <- log2_convert(ft_od)

# Normalize to treatment
#ft_log_norm <- normalize_table_to_treatment(metabolites_table = ft_log, metadata_table = LC_metadata, samples_norm = colnames(select(ft_log, contains("PBS"))))
ft_od_log_tnorm <- normalize_table_to_treatment(feature_table = ft_od_log, metadata_table = LC_metadata,
                                            grouping_variable = "syncom", samples_group_to_norm = "PBS")

write.csv(x = ft_od_log_tnorm, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/ft_od_log_tnorm.csv")

graph_metabolites(ft_od_log_tnorm, y1 = -10, y2 = 11, dotsize = 25, binwidth = 0.01)

##### log transformed + treatment-based log fold change

# Log transform 
ft_log <- log2_convert(ft_subset)

# Normalize to treatment
#ft_log_norm <- normalize_table_to_treatment(metabolites_table = ft_log, metadata_table = LC_metadata, samples_norm = colnames(select(ft_log, contains("PBS"))))
ft_log_tnorm <- normalize_table_to_treatment(feature_table = ft_log, metadata_table = LC_metadata,
                                                grouping_variable = "syncom", samples_group_to_norm = "PBS")

write.csv(x = ft_log_tnorm, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/ft_log_tnorm.csv")
# Graph
graph_metabolites(ft_log_tnorm, y1 = -10, y2 = 5, dotsize = 25, binwidth = 0.01)

  
# DO PCA
scaled_ft <- scale_0_1(fia_table, scale_by = "variable")

write.csv(x = ft_od_log_tnorm, file = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/ft_od_log_tnorm.csv")

scaled_ft_aa <- scale_0_1(ft_subset, scale_by = "variable")

ft_pca_2(scaled_ft, metadata_table = LC_metadata, grouping_col = "syncom")

ft_pca_2(scaled_ft_aa, metadata_table = LC_metadata, grouping_col = "syncom")
