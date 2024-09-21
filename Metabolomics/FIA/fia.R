source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")

# Read tables
fia_pos_table <- read_fia_table(table_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_SynComBatch1_2/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx")

fia_metadata_df <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_SynComBatch1_2/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx", sheet = "metadata", col_names = TRUE)

# Normalize by OD of SynComs.
fia_norm <- normalize_by_od(feature_table = fia_pos_table, metadata_table = fia_metadata_df, od_column = "OD")

# Remove highly variable metabolites.
filtered_fia <- filter_by_error(feature_table = fia_norm, metadata_table = fia_metadata_df, grouping_var = "SynCom", error_threshold = 50)

# Do PCA
fia_pca(feature_table = filtered_fia, metadata_table = fia_metadata_df, grouping_col = "SynCom")
