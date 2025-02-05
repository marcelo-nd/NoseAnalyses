library(dplyr)

source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

# Example dataframe
df <- data.frame(
  Species = c("Species_A", "Species_B", "Species_C"),
  sample1_r1 = c(5, 10, 15),
  sample1_r2 = c(6, 9, 14),
  sample1_r3 = c(4, 11, 16),
  sample2_r1 = c(20, 15, 10),
  sample2_r2 = c(22, 13, 8),
  sample2_r3 = c(19, 16, 12),
  sample3_r1 = c(8, 12, 20),
  sample3_r2 = c(9, 11, 21),
  sample3_r3 = c(7, 13, 19)
)

# Compute mean abundance per species for each sample
df_means <- df %>%
  # Pivot to long format for easier manipulation
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
  # Extract sample IDs without replicate information
  mutate(Sample_ID = sub("_r[0-9]+$", "", Sample)) %>%
  # Group by species and sample to calculate mean abundance
  group_by(Species, Sample_ID) %>%
  summarize(Mean_Abundance = mean(Abundance), .groups = "drop") %>%
  # Pivot back to wide format
  pivot_wider(names_from = Sample_ID, values_from = Mean_Abundance)

# View the resulting dataframe
print(df_means)

df_means <- as.data.frame(df_means)

print(df_means)

row.names(df_means) <- df_means$Species

df_means <- df_means[2:length(df_means)]

print(df_means)

barplot_from_feature_table(df_means, legend_cols = 1)

###################################
# Convert data to relative abundance

df_rel <- df

df_rel[, -1] <- df_rel[, -1] / colSums(df_rel[, -1])

### relative abundance averages

df_means_rel <- df_rel %>%
  # Pivot to long format for easier manipulation
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
  # Extract sample IDs without replicate information
  mutate(Sample_ID = sub("_r[0-9]+$", "", Sample)) %>%
  # Group by species and sample to calculate mean abundance
  group_by(Species, Sample_ID) %>%
  summarize(Mean_Abundance = mean(Abundance), .groups = "drop") %>%
  # Pivot back to wide format
  pivot_wider(names_from = Sample_ID, values_from = Mean_Abundance)

#### Relative abundances taking into account the other species

cumulative_subtraction <- function(df) {
  # Get numeric part of the dataframe (excluding species names)
  abundance_matrix <- as.matrix(df[, -1])
  
  # Apply cumulative subtraction for each column (sample)
  for (col in 1:ncol(abundance_matrix)) {
    remaining <- 1 # Start with 1 for each sample
    for (row in 1:nrow(abundance_matrix)) {
      current_value <- abundance_matrix[row, col]
      abundance_matrix[row, col] <- remaining - current_value
      remaining <- abundance_matrix[row, col] # Update remaining for next iteration
    }
  }
  
  # Replace values back into the dataframe
  df[, -1] <- abundance_matrix
  return(df)
}

# Apply the function
adjusted_df <- cumulative_subtraction(df_means_rel)


df_expl <- data.frame(
  Species = c("Species_A", "Species_B", "Species_C"),
  Sample1 = c(0.1667, 0.3333, 0.5000),
  Sample2 = c(0.1667, 0.3333, 0.5000),
  Sample3 = c(0.2000, 0.3000, 0.5000)
)


cumulative_subtraction(df_expl)
########################################
df <- data.frame(
  Species = c("Species_A", "Species_B", "Species_C"),
  Sample1 = c(0.1667, 0.3333, 0.5000),
  Sample2 = c(0.1667, 0.3333, 0.5000),
  Sample3 = c(0.2000, 0.3000, 0.5000)
)

cumulative_half_addition <- function(df) {
  # Get numeric part of the dataframe (excluding species names)
  abundance_matrix <- as.matrix(df[, -1])
  
  # Apply calculations for each column (sample)
  for (col in 1:ncol(abundance_matrix)) {
    for (row in nrow(abundance_matrix):1) { # Iterate from bottom to top
      if (row == nrow(abundance_matrix)) {
        # Last row: itself divided by 2
        abundance_matrix[row, col] <- abundance_matrix[row, col] / 2
      } else {
        # Other rows: cumulative sum of rows below + itself/2
        abundance_matrix[row, col] <- sum(abundance_matrix[(row + 1):nrow(abundance_matrix), col]) +
          (abundance_matrix[row, col] / 2)
      }
    }
  }
  
  # Replace values back into the dataframe
  df[, -1] <- abundance_matrix
  return(df)
}

cumulative_half_addition(df)


##### Calculate relative abundance
df_rel <- df

df_rel[, -1] <- df_rel[, -1] / colSums(df_rel[, -1])

##### Calculate cumulative half addition
df_hcs <- cumulative_half_addition(df_rel)

##### Create long format dataframe
df_long <- df_hcs %>%
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance")

df_long <- df_rel %>%
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Sample_ID = sub("_r[0-9]+$", "", Sample),  # Extract sample name without replicate info
    Replicate = sub(".*_r", "", Sample)        # Extract replicate number
  )

##### Plot
ggplot2::ggplot(df_long, ggplot2::aes(x=Sample, y=Abundance, fill=Species, color = Species)) + 
  ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
  geom_point(ggplot2::aes(shape = Species), size = 4)

ggplot2::ggplot(df_long, ggplot2::aes(x=Sample, y=Abundance, fill=Species, shape = Species)) + 
  geom_point(size = 4)

#### Dot barplot
# Transform to long format for easier plotting
df_long <- df_rel %>%
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Sample_ID = sub("_r[0-9]+$", "", Sample),  # Extract sample name without replicate info
    Replicate = sub(".*_r", "", Sample)        # Extract replicate number
  )

# Calculate relative abundance per sample
df_long <- df_long %>%
  group_by(Sample_ID) %>%
  mutate(Rel_Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# Create dot plot
ggplot(df_long, aes(x = Sample_ID, y = Abundance, color = Species, shape = Replicate)) +
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +  # Add dots with jitter
  labs(x = "Sample", y = "Relative Abundance", color = "Species", shape = "Replicate") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



ggplot2::ggplot(df_long, ggplot2::aes(x=Sample_ID, y=Abundance, fill=Species, color = Species)) + 
  ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
  geom_point(size = 4)

##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################

df <- data.frame(
  Species = c("Species_A", "Species_B", "Species_C"),
  Sample1 = c(0.1667, 0.3333, 0.5000),
  Sample2 = c(0.1667, 0.3333, 0.5000),
  Sample3 = c(0.2000, 0.3000, 0.5000)
)

cumulative_half_addition <- function(df) {
  # Get numeric part of the dataframe (excluding species names)
  abundance_matrix <- as.matrix(df[, -1])
  
  df2 <- as.matrix(df[, -1])
  
  # Apply calculations for each column (sample)
  for (col in 1:ncol(abundance_matrix)) {
    for (row in nrow(abundance_matrix):1) { # Iterate from bottom to top
      if (row == nrow(abundance_matrix)) {
        # Last row: itself divided by 2
        if (df2[row, col] == 0) {
          abundance_matrix[row, col] <- 0
        } else{
          abundance_matrix[row, col] <- df2[row, col] / 2
        }
      } else { # Other rows: cumulative sum of rows below + itself/2
        if (df2[row, col] == 0) {
          abundance_matrix[row, col] <- 0
        } else {
          abundance_matrix[row, col] <- sum(df2[(row + 1):nrow(df2), col]) +
            (df2[row, col] / 2)
        }
      }
    }
  }
  
  # Replace values back into the dataframe
  df[, -1] <- abundance_matrix
  return(df)
}

##### Calculate relative abundance
df_rel <- df
df_rel <- sc4
df_rel <- rownames_to_column(df_rel, var = "Species")

calculate_relative_abundance <- function(df) {
  # Ensure the first column is species names
  species <- df[, 1]
  
  # Extract the numeric sample data
  sample_data <- df[, -1]
  
  # Calculate relative abundance
  relative_abundance <- sweep(sample_data, 2, colSums(sample_data), "/")
  
  # Combine species names back with the relative abundance data
  result <- cbind(Species = species, relative_abundance)
  
  # Return the result as a dataframe
  return(as.data.frame(result))
}

df_rel <- calculate_relative_abundance(df_rel)

##### Calculate cumulative half addition
df_hcs <- cumulative_half_addition(df_rel)

##### Create long format dataframe
# For dots
df_long_dots <- df_hcs %>%
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance")

ggplot2::ggplot(df_long_dots, ggplot2::aes(x=Sample, y=Abundance, fill=Species, color = Species)) + 
  geom_point(size = 4)

# For bars
df_long_bars <- df_rel %>%
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance")

ggplot2::ggplot(df_long_bars, ggplot2::aes(x=Sample, y=Abundance, fill=Species, color = Species)) + 
  ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE)

df_long_joined <- df_long_bars

df_long_joined$Abundance_dots <- df_long_dots$Abundance

ggplot2::ggplot(df_long_joined) + 
  ggplot2::geom_bar(aes(x=Sample, y=Abundance, fill=Species, color = Species), position="fill", stat="identity", show.legend = TRUE) +
  geom_point(aes(x=Sample, y=Abundance_dots, shape = Species), size = 4)


### Now lets put together replicate data
### Calculate the means of each sample
df_means_rel <- df_rel %>%
  # Pivot to long format for easier manipulation
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
  # Extract sample IDs without replicate information
  mutate(Sample_ID = sub("_R[0-9]+$", "", Sample)) %>%
  # Group by species and sample to calculate mean abundance
  group_by(Species, Sample_ID) %>%
  summarize(Mean_Abundance = mean(Abundance), .groups = "drop") #%>%
  # Pivot back to wide format
  #pivot_wider(names_from = Sample_ID, values_from = Mean_Abundance)

##### Plot
ggplot2::ggplot(df_means_rel, ggplot2::aes(x=Sample_ID, y=Mean_Abundance, fill=Species, color = Species)) + 
  ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE)

### Put together points

#df_long <- df_rel %>%
#  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
#  mutate(
#    Sample_ID = sub("_r[0-9]+$", "", Sample),  # Extract sample name without replicate info
#    Replicate = sub(".*_r", "", Sample)        # Extract replicate number
#  )

df_long_dots2 <- df_long_dots %>%
  mutate(
    Sample_ID = sub("_R[0-9]+$", "", Sample),  # Extract sample name without replicate info
    Replicate = sub(".*_R", "", Sample)        # Extract replicate number
  ) %>%
  filter(Abundance != 0)

ggplot2::ggplot(df_long_dots2, ggplot2::aes(x=Sample_ID, y=Abundance, fill=Species, colour = Species)) + 
  geom_point(size = 4)


(plot2 <- ggplot(NULL) + 
    ggplot2::geom_bar(data = df_means_rel, aes(x=Sample_ID, y=Mean_Abundance, fill=Species, color = Species), position="fill", stat="identity", show.legend = TRUE) +
    geom_point(data = df_long_dots2, aes(x=Sample_ID, y=Abundance, shape = Species), size = 3, color = "white", stroke = 1.3) +
    scale_shape_manual(values = c(0, 1, 17, 3, 4, 5, 6, 20, 8, 9, 10)) +
    ggthemes::theme_tufte()
)

barplot_with_replicates(sc4)

barplot_with_replicates(sc7)

barplot_with_replicates(sc9)

barplot_with_replicates(sc11)

barplot_with_replicates(sc13)

barplot_with_replicates(sc20)

barplot_with_replicates(sc24)

barplot_with_replicates(sc25)

barplot_with_replicates(sc32)

barplot_with_replicates(sc43)

barplot_with_replicates(sc53)






