library(ggplot2)
library(ggfortify)
library(plotly)

# Read fia positive mode
fia_pos_df <- readxl::read_excel(path = "C:/Users/marce/Desktop/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx", sheet = "pos", col_names = TRUE)

# Read fia negative mode
fia_pos_df <- readxl::read_excel(path = "C:/Users/marce/Desktop/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx", sheet = "neg", col_names = TRUE)

# Read metadata
fia_metadata_df <- readxl::read_excel(path = "C:/Users/marce/Desktop/20240523_FIA_Fc_SynComBatch1_2_fcexport.xlsx", sheet = "metadata", col_names = TRUE)


# Retain only columns necessary for PCA

fia_pos_df_2 <- cbind(fia_pos_df[, 2], fia_pos_df[, 5:ncol(fia_pos_df)])

# Transpose table

fia_pos_df_2 <- t(fia_pos_df_2)

# Set column names as metabolites names
colnames(fia_pos_df_2) <- fia_pos_df_2[1,]

# Remove 
fia_pos_df_2 <- fia_pos_df_2[2:nrow(fia_pos_df_2),]

# Check number of rows is the same as number of samples and that they are in the same order in Metadata and feature table.
nrow(fia_pos_df_2)

rownames(fia_pos_df_2) == fia_metadata_df$Sample

# Reconvert to dataframe
fia_pos_df_2 <- as.data.frame(fia_pos_df_2)

# Transform to numeric all columns
fia_pos_df_2[,1:ncol(fia_pos_df_2)] <- sapply(fia_pos_df_2[,1:ncol(fia_pos_df_2)],as.numeric)

# Check the data types of columns
sapply(fia_pos_df_2, class)

# Add metadata to same dataframe

joined_df <- cbind(fia_metadata_df, fia_pos_df_2)

colnames(joined_df) <- make.names(colnames(joined_df), unique=TRUE)

# Normalize by OD of SynComs

normalize_by_od <- function(dataframe, select_by = "range", select_info){
  if (select_by == "range" || select_by == "names") {
    df_normalized <- dataframe %>%
      mutate(across(all_of(select_info), ~ . / OD)) %>% # to divide based on names list or indices
      #mutate(across(starts_with("var"), ~ . / div_value)) %>% # to divide variables by pattern in name
      #select(-div_value)  # Remove the division value column if not needed
    # Return the normalized dataframe
    return(df_normalized)
  }
  else if (select_by == "numeric") {
    df_normalized <- dataframe %>%
      mutate(across(where(is.numeric), ~ . / div_value)) %>% # to divide all numeric variables
      # Return the normalized dataframe
      return(df_normalized)
  }
  else if (select_by == "pattern") {
    df_normalized <- dataframe %>%
      mutate(across(starts_with(select_info), ~ . / div_value)) %>% # to divide variables by pattern in name e.g. "var"
      # Return the normalized dataframe
      return(df_normalized)
  }
  else{
    print("select_by value not valid")
  }
}

fia_norm <- normalize_by_od(joined_df, select_info = c(5: ncol(joined_df)))

print(head(fia_norm))

# Remove highly variable metabolites

filter_by_error <- function(dataframe, error_threshold = 50){
  
  # Step 1: Calculate the error for each variable per type
  errors <- dataframe %>%
    group_by(SynCom) %>%
    summarise(across(where(is.numeric), ~ (sd(.) / mean(.)) * 100, .names = "error_{col}"))
  
  # View the errors dataframe
  #print(head(errors))
  
  # Step 2: Average the error for each variable for all types
  avg_errors <- errors %>%
    summarise(across(starts_with("error_"), mean, na.rm = TRUE))
  
  # View the average errors dataframe
  #print(head(avg_errors))
  
  variables_to_keep <- names(avg_errors)[avg_errors <= threshold]
  variables_to_keep <- gsub("error_", "", variables_to_keep)
  
  # Keep the metadata columns and the variables with error below the threshold
  df_filtered <- dataframe %>%
    select(all_of(c("Sample", "SynCom", "Time", "OD", variables_to_keep)))
  
  # View the filtered dataframe
  return(df_filtered)
}

filtered_fia <- filter_by_error(fia_norm, error_threshold = 50)

# Do PCA
fia_pos_pca <- prcomp(fia_pos_df_2, scale. = TRUE)

fia_pos_pca <- prcomp(filtered_fia[, 5:ncol(filtered_fia)], scale. = TRUE)

# Plot PCA with samples coloured by SynCom
p2 <- ggplot2::autoplot(fia_pos_pca, data = na.omit(fia_metadata_df), colour = 'SynCom') +
  #Size and alpha just for fun
  scale_color_manual(values = get_palette(60))

p2
