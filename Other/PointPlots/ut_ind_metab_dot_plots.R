library(dplyr)
library(ggplot2)

TURKEY.TEST.ALL.SIGNIFIANT.FEATURES <- read.csv("C:/Users/marce/OneDrive - UT Cloud/PointPlots/TURKEY TEST ALL SIGNIFIANT FEATURES.csv", header=FALSE, sep=";")

# Unique groups
unique(TURKEY.TEST.ALL.SIGNIFIANT.FEATURES$V2)

X2023_09_26CLEANED_DATA_WITH_MD <- readr::read_csv("C:/Users/marce/OneDrive - UT Cloud/PointPlots/2023-09-26CLEANED_DATA_WITH_MD.csv")

# Extract data for two groups
extract_features_comparison <- function(comparison_name, clean_data, tukey_test_table){
  sig_features <- TURKEY.TEST.ALL.SIGNIFIANT.FEATURES[TURKEY.TEST.ALL.SIGNIFIANT.FEATURES$V2 == comparison_name,]$V1
  #asulam_features <- asulam$V1
  temp_df <- dplyr::select(clean_data, any_of(c("Row.names", "ATTRIBUTE_GROUPS", sig_features)))
  temp_df <- temp_df[2:length(temp_df)]
  temp_df_t <- as.data.frame(t(temp_df))
  
  colnames(temp_df_t) <- temp_df$ATTRIBUTE_GROUPS
  
  colnames(temp_df_t) <- make.unique(colnames(temp_df_t), sep = "_")
  
  temp_df_t <- temp_df_t[2:nrow(temp_df_t),]
  
  # Convert to numeric
  i = colnames(temp_df_t)
  temp_df_t[ , i] <- apply(temp_df_t[ , i], 2,  # Specify own function within apply
                           function(x) as.numeric(as.character(x)))
  
  temp_df_t$Metabolite <- row.names(temp_df_t)
  return(temp_df_t)
}

#Asulam

asulam_features <- extract_features_comparison(comparison_name = "Asulam_24-Asulam_1",
                                      clean_data = X2023_09_26CLEANED_DATA_WITH_MD,
                                      tukey_test_table = TURKEY.TEST.ALL.SIGNIFIANT.FEATURES)

asulam_values <- dplyr::select(asulam_features, starts_with("Asulam"), starts_with("Metabolite"))

asulam_metabolites_norm <- normalize_table_to_treatment(asulam_values, "Asulam_1")

asu1 <- gather_samples(asulam_metabolites_norm, sample.names = c("Asulam_1", "Asulam_1_1", "Asulam_1_2"), sample_type = "Asulam_day1")

asu24 <- gather_samples(asulam_metabolites_norm, sample.names = c("Asulam_24", "Asulam_24_1", "Asulam_24_2"), sample_type = "Asulam_day24")

asu_data <- rbind(asu1, asu24)

graph_metabolites(asu_data, y1 = -5, y2 = 5, dotsize = 15, binwidth = 0.01, ylab = "Norm. intensity")

# Glyphosate

Glyphosate_features <- extract_features_comparison(comparison_name = "Glyphosate_24-Glyphosate_1",
                                               clean_data = X2023_09_26CLEANED_DATA_WITH_MD,
                                               tukey_test_table = TURKEY.TEST.ALL.SIGNIFIANT.FEATURES)

Glyphosate_values <- dplyr::select(Glyphosate_features, starts_with("Glyphosate"), starts_with("Metabolite"))

Glyphosate_metabolites_norm <- normalize_table_to_treatment(Glyphosate_values, "Glyphosate_1")

gly1 <- gather_samples(Glyphosate_metabolites_norm, sample.names = c("Glyphosate_1", "Glyphosate_1_1", "Glyphosate_1_2"), sample_type = "Glyphosate_day1")

gly24 <- gather_samples(Glyphosate_metabolites_norm, sample.names = c("Glyphosate_24", "Glyphosate_24_1", "Glyphosate_24_2"), sample_type = "Glyphosate_day24")

gly_data <- rbind(gly1, gly24)

graph_metabolites(gly_data, y1 = -5, y2 = 5, dotsize = 15, binwidth = 0.01, ylab = "Norm. intensity")


# Fosfomycin

Fosfomycin_features <- extract_features_comparison(comparison_name = "Fosfomycin_24-Fosfomycin_1",
                                                   clean_data = X2023_09_26CLEANED_DATA_WITH_MD,
                                                   tukey_test_table = TURKEY.TEST.ALL.SIGNIFIANT.FEATURES)

Fosfomycin_values <- dplyr::select(Fosfomycin_features, starts_with("Fosfomycin"), starts_with("Metabolite"))

Fosfomycin_metabolites_norm <- normalize_table_to_treatment(Fosfomycin_values, "Fosfomycin_1")

fos1 <- gather_samples(Fosfomycin_metabolites_norm, sample.names = c("Fosfomycin_1", "Fosfomycin_1_1", "Fosfomycin_1_2"), sample_type = "Fosfomycin_day1")

fos24 <- gather_samples(Fosfomycin_metabolites_norm, sample.names = c("Fosfomycin_24", "Fosfomycin_24_1", "Fosfomycin_24_2"), sample_type = "Fosfomycin_day24")

fos_data <- rbind(fos1, fos24)

fos_graph <- graph_metabolites(fos_data, y1 = -5, y2 = 5, dotsize = 15, binwidth = 0.01, ylab = "Norm. intensity")

fos_graph + theme(axis.text=element_text(size=5))


# Imidazole
Imidazole_features <- extract_features_comparison(comparison_name = "Imidazole_24-Imidazole_1",
                                                   clean_data = X2023_09_26CLEANED_DATA_WITH_MD,
                                                   tukey_test_table = TURKEY.TEST.ALL.SIGNIFIANT.FEATURES)

Imidazole_values <- dplyr::select(Imidazole_features, starts_with("Imidazole"), starts_with("Metabolite"))

Imidazole_metabolites_norm <- normalize_table_to_treatment(Imidazole_values, "Imidazole_1")

imi1 <- gather_samples(Imidazole_metabolites_norm, sample.names = c("Imidazole_1", "Imidazole_1_1", "Imidazole_1_2"), sample_type = "Imidazole_day1")

imi24 <- gather_samples(Imidazole_metabolites_norm, sample.names = c("Imidazole_24", "Imidazole_24_1", "Imidazole_24_2"), sample_type = "Imidazole_day24")

imi_data <- rbind(imi1, imi24)

imi_graph <- graph_metabolites(imi_data, y1 = -5, y2 = 5, dotsize = 15, binwidth = 0.01, ylab = "Norm. intensity")

imi_graph + theme(axis.text=element_text(size=5))

