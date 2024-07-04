Sys.setenv(LANG = "en")

library(ggplot2)
library(ggfortify)
library(plotly)
library(dplyr)

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

library(PCAtools)

source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")

# Read fia positive mode
fia_df <- readxl::read_excel(path = "C:/Users/marce/Desktop/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx", sheet = "pos", col_names = TRUE)

# Read fia negative mode
fia_df <- readxl::read_excel(path = "C:/Users/marce/Desktop/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx", sheet = "neg", col_names = TRUE)

# Read metadata
fia_metadata_df <- readxl::read_excel(path = "C:/Users/marce/Desktop/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx", sheet = "metadata", col_names = TRUE)

# Retain only columns necessary for PCA

fia_df_2 <- cbind(fia_df[, 2], fia_df[, 5:ncol(fia_df)])

# Transpose table

fia_df_t <- t(fia_df_2)

# Set column names as metabolites names
colnames(fia_df_t) <- fia_df_t[1,]

# Remove 
fia_df_t <- fia_df_t[2:nrow(fia_df_t),]

# Check number of rows is the same as number of samples and that they are in the same order in Metadata and feature table.
nrow(fia_df_t)

rownames(fia_df_t) == fia_metadata_df$Sample

# Reconvert to dataframe
fia_df_t <- as.data.frame(fia_df_t)

# Transform to numeric all columns
fia_df_t[,1:ncol(fia_df_t)] <- sapply(fia_df_t[,1:ncol(fia_df_t)],as.numeric)

# Check the data types of columns
head(sapply(fia_df_t, class))

# Add metadata to same dataframe

joined_df <- cbind(fia_metadata_df, fia_df_t)

colnames(joined_df) <- make.names(colnames(joined_df), unique=TRUE)

# Normalize by OD of SynComs

fia_norm <- normalize_by_od(joined_df, select_info = c(5: ncol(joined_df)))

print(head(fia_norm))

# Remove highly variable metabolites

filtered_fia <- filter_by_error(fia_norm, error_threshold = 50)

fia_pos_pca <- prcomp(filtered_fia[, 5:ncol(filtered_fia)], scale. = TRUE)


fia_metadata_df <- as.data.frame(fia_metadata_df)

# Plot PCA with samples coloured by SynCom
p2 <- ggplot2::autoplot(fia_pos_pca, data = na.omit(fia_metadata_df), colour = 'SynCom') +
  # change color scale
  scale_color_manual(values = get_palette(60))

p2



###############################################################################
# Reprocess Fia df

fia_df_3 <- t(filtered_fia[,5:ncol(filtered_fia)])

# Check number of rows is the same as number of samples and that they are in the same order in Metadata and feature table.
ncol(fia_df_3)

colnames(fia_df_3) == fia_metadata_df$Sample

# Metadata
fia_metadata_df2 <- fia_metadata_df
fia_metadata_df2 <- as.data.frame(fia_metadata_df2)

rownames(fia_metadata_df2) <- fia_metadata_df2$Sample

fia_pos_pca <- PCAtools::pca(fia_df_3, scale = TRUE, metadata =fia_metadata_df2, transposed = FALSE)

p2 <- PCAtools::biplot(fia_pos_pca, showLoadings = TRUE, ntopLoadings = 2, lab = NULL, colby = "SynCom",
                       legendPosition = "right", legendLabSize = 9, legendIconSize = 2, pointSize = 2,
                       colkey = get_palette(nColors = 26))

p2
