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




