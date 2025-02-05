source("C:/Users/marce/Documents/GitHub/microbiome-help/functionalDataWrangling.R")

library(dplyr)

feature_table <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/annotated_QuantTable_scaled.csv")

feature_table <- feature_table[rev(order(feature_table[["X"]])), ]

#names(feature_table)[names(feature_table) == 'X'] <- 'filename'

row.names(feature_table) <- feature_table$X

feature_table <- feature_table[-1]


metadata <- read.csv("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/SC100_metadata_noqcs_nosinStrs.csv")

metadata <- metadata[7:nrow(metadata),]

metadata <- metadata[rev(order(metadata[["filename"]])), ]

names(metadata)[names(metadata) == 'filename'] <- 'Sample'

rownames(feature_table) == metadata$Sample


# Remove highly variable metabolites (intrareplicate variability).
filtered_ft <- filter_by_error(feature_table = feature_table, metadata_table = metadata, grouping_var = "ATTRIBUTE_Sample", error_threshold = 25)
rownames(filtered_ft) <- rownames(feature_table)

# Do PCA
ft_pca(feature_table = filtered_ft, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", encircle = FALSE)

ft_pca(feature_table = filtered_ft, metadata_table = metadata, grouping_col = "ATTRIBUTE_Sample", encircle = TRUE)



# PCA with vegan
library(vegan)
library(ggplot2)
library(dplyr)


# Compute distance matrix (Bray-Curtis, but you can change to another method)
dist_matrix <- vegdist(filtered_ft, method = "euclidean")

# Perform PCoA
pcoa_results <- cmdscale(dist_matrix, k = 2, eig = TRUE)

pcoa_df <- as.data.frame(pcoa_results$points)
colnames(pcoa_df) <- c("PC1", "PC2")  # Rename axes
pcoa_df$Sample <- rownames(feature_table)  # Add sample names

print(pcoa_df$Sample == metadata$Sample)

# Merge with metadata to get Sample_type
pcoa_df <- left_join(pcoa_df, metadata, by = c("Sample" = "Sample"))

# Plot PCoA
colour_palette <- get_palette(nColors = 50)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = ATTRIBUTE_Sample, , shape = ATTRIBUTE_Time)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values=colour_palette) +
  labs(title = "PCoA Plot (Bray-Curtis)",
       x = "PCoA 1",
       y = "PCoA 2",
       color = "Sample Type") +
  theme(legend.position = "right")
