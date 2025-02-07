source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")

feature_table <- read_ft_1("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/annotated_QuantTable_scaled.csv",
                             sort_by_names = TRUE)

#feature_table <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/annotated_QuantTable_scaled.csv")

#feature_table <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_scaled_an.csv")

#feature_table <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_scaled.csv")

#feature_table <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_imp_an.csv")

#feature_table <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_imp.csv")

#feature_table <- feature_table[rev(order(feature_table[["X"]])), ]

#feature_table <- feature_table[order(feature_table[["X"]]), ]

#names(feature_table)[names(feature_table) == 'X'] <- 'filename'

#row.names(feature_table) <- feature_table$filename

#feature_table <- feature_table[-1]

metadata <- read_metadata("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/SC100_metadata_noqcs_nosinStrs.csv",
                            sort_by_names = TRUE)

metadata <- metadata[7:nrow(metadata),]

# Remove highly variable metabolites (intrareplicate variability).
filtered_ft <- filter_by_error(feature_table = feature_table, metadata_table = metadata, grouping_var = "ATTRIBUTE_Sample", error_threshold = 25)
rownames(filtered_ft) <- rownames(feature_table)

# Do PCA
ft_pca(feature_table = feature_table, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", encircle = FALSE)

ft_pca(feature_table = filtered_ft, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", encircle = TRUE)

# PCA with vegan
ft_pca_2(feature_table = feature_table, metadata_table = metadata, dist_method = "euclidean",
         grouping_col = "ATTRIBUTE_Sample", p_shape = "ATTRIBUTE_Time")


