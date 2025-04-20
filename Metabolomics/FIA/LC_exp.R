source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")


# Read data
fia_table <- read_fia_table("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_SC_exp/20250115_FIA-Data_lc_SCs_exp.xlsx",
                            sort_table = TRUE, fix_names = FALSE)

# Read targeted analyses table
ft <- read_targeted_table("D:/2_OtherProjects/LauraSynComs/20250408_ratioList_final.xlsx", sheet = "Ratios",
                     sort_table = TRUE, in_col = 9)

# Read Metadata
LC_metadata <- read_metadata_xls(path = "D:/2_OtherProjects/LauraSynComs/Metabolomics tubes annotation - Nov 2024_2.xlsx",
                                 sheet = "Tabelle1", sort_table = TRUE)

write.csv(LC_metadata, "D:/2_OtherProjects/LauraSynComs/metadata_table.csv")

# Imputing data table
ft_imp <- ft_imputation(ft)

# subset only aminoacids
# ft_aas <- extract_features_comparison(feature_table = ft_imp, sig_features = c("L-Alanine[M+H]+", "L-Arginine[M+H]+", "L-Asparagine[M+H]+",
#                                                                                   "L-Aspartic acid[M+H]+", "L-Cystine[M+H]+", "L-Glutamic acid[M+H]+",
#                                                                                   "L-Glutamine[M+H]+", "L-Histidine[M+H]+", "L-Isoleucine[M+H]+",
#                                                                                   "L-Leucine[M+H]+", "L-Lysine[M+H]+", "L-Methionine[M+H]+",
#                                                                                   "L-Phenylalanine[M+H]+", "L-Proline[M+H]+", "L-Serine[M+H]+",
#                                                                                   "L-Threonine[M+H]+", "L-Tryptophan[M+H]+", "L-Tyrosine[M+H]+",
#                                                                                   "L-Valine[M+H]+", "Glycine[M+H]+"), columns_to_preserve = NULL)

ft_aas <- extract_features_comparison(feature_table = ft_imp, sig_features = c("Alanine", "Arginine", "Asparagine", "Glutamic acid", "Glutamine",
                                                                                "Glycine", "Isoleucine", "Leucine", "Lysine", "Methionine",
                                                                                "Phenylalanine", "Proline", "Serine", "Threonine/Homoserine",
                                                                                "Tryptophan", "Tyrosine", "Valine"), columns_to_preserve = NULL)

# Read LCMS table
lcms_table <- read_lcms_table("D:/2_OtherProjects/LauraSynComs/20250408_ratioList_final.xlsx")

write.csv(lcms_table, "D:/2_OtherProjects/LauraSynComs/lcms_table.csv")

# Chose table to use
ft <- ft_imp
ft <- ft_aas
ft <- ft_aas_scaled
ft <- lcms_table

# Scale feature table for PCA
ft_scaled <- scale_0_1(ft, scale_by = "variable") # What is variable?

# Do PCAs
color_p = c("#89617E", "#7AAA83", "#AFABAB", "#F4D171", "#3A3937", "#17678B")

# PCA
ft_pca <- ft_pca_1(ft_scaled, metadata_table = LC_metadata, grouping_col = "syncom", color_palette = color_p, encircle = TRUE)


ft_pca

###### Fold change plots
# Normalize to OD
ft_od <- normalize_by_od(feature_table = ft, metadata_table = LC_metadata, grouping_variable = "syncom", samples_group_to_exclude = "SNM")

# log fold calculation, no OD norm
ft_lf <- log2_fold(ft, metadata_table = LC_metadata,
                     grouping_variable = "syncom", samples_group_to_norm = "SNM")

# log fold calculation, OD norm
ft_lf_od <- log2_fold(ft_od, metadata_table = LC_metadata,
                   grouping_variable = "syncom", samples_group_to_norm = "SNM")

# Fold change plot, no OD norm
graph_metabolites(ft_lf, y1 = -10, y2 = 10, dotsize = 25, binwidth = 0.01)

# Fold change plot, OD norm
graph_metabolites(ft_lf_od, y1 = -10, y2 = 10, dotsize = 25, binwidth = 0.01)
