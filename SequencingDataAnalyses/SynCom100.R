# Load helping functions
source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")

# Screening

otu_table_screening <- read.csv("F:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1)

colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8",
                                   "SC9", "SC10", "SC11", "SC12", "SC13", "SC14", "SC15",
                                   "SC16", "SC17", "SC18", "SC19", "SC20", "SC21", "SC22",
                                   "SC23", "SC24", "SC25", "SC26", "SC27", "SC28", "SC29",
                                   "SC30", "SC31", "SC32", "SC33", "SC34", "SC35", "SC36",
                                   "SC37", "SC38", "SC39", "SC40", "SC41", "SC42", "SC43",
                                   "SC44", "SC45", "SC46", "SC47", "SC48", "SC49", "SC50",
                                   "SC51", "SC52", "SC53")

otu_table_screening_filt <- filter_otus_by_counts_col_counts(otu_table_screening,
                                                             min_count = 10,
                                                             col_number = 1)

df_filtered <- otu_table_screening[, -39]

df_filtered <- df_filtered[, -34]

df_filtered <- df_filtered[, -5]

df_filtered <- df_filtered[-3,]

df_filtered <- df_filtered[-20,]

df_filtered <- df_filtered[apply(df_filtered, 1, function(row) sum(row >= 10) >= 1), ]

otu_table_adjusted <- df_filtered
otu_table_adjusted[6, 25:29] <- otu_table_adjusted[7, 25:29] / 5

barplot_from_feature_table(otu_table_adjusted)


##### Timepoints

otu_table_SC_tp <- read.csv("F:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv", row.names=1)

cf <- colnames(otu_table_SC_tp)

ordered_names <- cf[order(nchar(cf), cf)]

otu_table_SC_tp <- otu_table_SC_tp[, ordered_names]

colnames(otu_table_SC_tp) <- c("SC4_TF_R1", "SC4_TF_R2", "SC4_TF_R3",
                              "SC7_TF_R1", "SC7_TF_R2", "SC7_TF_R3",
                              "SC13_TF_R1", "SC13_TF_R2", "SC13_TF_R3",
                              "SC20_TF_R1", "SC20_TF_R2", "SC20_TF_R3",
                              "SC24_TF_R1", "SC24_TF_R2", "SC24_TF_R3",
                              "SC25_TF_R1", "SC25_TF_R2", "SC25_TF_R3",
                              "SC32_TF_R1", "SC32_TF_R2", "SC32_TF_R3",
                              "SC43_TF_R1", "SC43_TF_R2", "SC43_TF_R3",
                              "SC47_TF_R1", "SC47_TF_R2", "SC47_TF_R3",
                              "SC53_TF_R1", "SC53_TF_R2", "SC53_TF_R3")

otu_table_SC_tp_filt <- filter_otus_by_counts_col_counts(otu_table_SC24,
                                                             min_count = 50,
                                                             col_number = 1)

otu_table_SC_tp <- otu_table_SC_tp[-2,]

otu_table_SC_tp <- otu_table_SC_tp[-9,]

barplot_from_feature_table(otu_table_SC_tp)


##### Ordered plots
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


######## First calculate the proportions of S. aureus.
df <- otu_table_adjusted %>% rownames_to_column(var = "Species")

#df <- df %>% rownames_to_column(var = "Species")

#df <- df[1:8,]

# Specify the species of interest
species_of_interest <- "Staphylococcus aureus"

# Calculate the total abundance per sample
total_abundance <- colSums(df[, -1])

# Filter the row of the species of interest and calculate its proportion with respect to total abundance
df_proportion <- df %>%
  filter(Species == species_of_interest) %>%
  select(-Species)

df_proportion <- df_proportion[1,]/total_abundance

ordered_samples <- df_proportion %>%
  unlist() %>%
  sort(decreasing = TRUE) %>%
  names()

df_long <- df %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")

df_long$Sample <- factor(df_long$Sample, levels = ordered_samples)

color_palette <- get_palette(nColors = 50)

otu_barplot <- ggplot2::ggplot(df_long, ggplot2::aes(x=Sample, y=Abundance, fill=Species)) + 
  ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::scale_fill_manual(values=colors_vec) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.title=ggplot2::element_text(size=10), 
                 legend.text=ggplot2::element_text(size=08))
otu_barplot


##### Plot with bars ordered by similarity
# Z-Scaling
df2 <- as.data.frame(scale(otu_table_adjusted))

# Min-max scaling
df2 <- otu_table_adjusted
normalize = function(x) (x- min(x))/(max(x) - min(x))
cols <- sapply(otu_table_adjusted, is.numeric)
df2[cols] <- lapply(otu_table_adjusted[cols], normalize)

df2 <- df2 %>% rownames_to_column(var = "Species")

df_t <- as.matrix(t(df2[, -1]))  # Exclude the "Species" column after moving it to row names

# Perform hierarchical clustering
d <- dist(df_t, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

# Get the order of samples based on clustering
ordered_samples_cluster <- colnames(df)[-1][hc$order]  # Remove "Species" again

# Update sample factor levels in the long-format data for ggplot
df_long$Sample <- factor(df_long$Sample, levels = ordered_samples_cluster)

colors_vec <- c("gold3", "olivedrab3", "blueviolet", "#CC79A7",
                "lightblue1","brown1", "#053f73", "darkorange3")

# Plot ordered by clustering similarity
otu_barplot <- ggplot2::ggplot(df_long, ggplot2::aes(x=Sample, y=Abundance, fill=Species)) + 
  ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::scale_fill_manual(values=colors_vec) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                 legend.title=ggplot2::element_text(size=10), 
                 legend.text=ggplot2::element_text(size=08))
otu_barplot
