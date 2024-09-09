library(dplyr)
library(ggplot2)

feature_table <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_Test_230824/20240822_FIA-Data_LC_Test_pos.xlsx")

feature_table <- read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/FIA/Results_LC_Test_230824/20240822_MergePos_LC_Test.xlsx")

LC_metadata <- read_excel("C:/Users/marce/OneDrive - UT Cloud/LC_metadata.xlsx")

# Extract data that matches a list of compounds froma  feature table.
extract_features_comparison <- function(comparison_name, feature_table, sig_features, columns_to_preserve = NULL){
  # Creat list of columns to extract
  if (!is.null(columns_to_preserve)) {
    columns_to_preserve <- c(columns_to_preserve, sig_features)
  }else{
    columns_to_preserve <- sig_features
  }
  
  print(columns_to_preserve)
  
  temp_df <- dplyr::select(feature_table, any_of(columns_to_preserve))
  print(head(temp_df))
  #temp_df <- temp_df[2:length(temp_df)]
  temp_df_t <- as.data.frame(t(temp_df))
  
  #colnames(temp_df_t) <- temp_df$ATTRIBUTE_GROUPS
  
  colnames(temp_df_t) <- make.unique(colnames(temp_df_t), sep = "_")
  
  #temp_df_t <- temp_df_t[2:nrow(temp_df_t),]
  
  # Convert to numeric
  i = colnames(temp_df_t)
  temp_df_t[ , i] <- apply(temp_df_t[ , i], 2,  # Specify own function within apply
                           function(x) as.numeric(as.character(x)))
  
  temp_df_t$Metabolite <- row.names(temp_df_t)
  return(temp_df_t)
}

log2_convert <- function(metabolites_table){
  metabolites_table_log2 <- metabolites_table %>% dplyr::select(-Metabolite) %>% log2()
  #metabolites_table_log2 <- metabolites_table %>% log2()
  #rownames(metabolites_table_log2) <- metabolites_table$Metabolite
  metabolites_table_log2 <- tibble::add_column(metabolites_table_log2, "Metabolite" = metabolites_table$Metabolite, .before = 1)
  return(metabolites_table_log2)
}

normalize_table_to_treatment <- function(metabolites_table, samples_norm){
  #metabolites_table_norm <- metabolites_table %>% dplyr::select(-Metabolite) %>% t()
  metabolites_table_norm <- metabolites_table %>% dplyr::select(-Metabolite)
  
  # Get the means of the treatment used to normalize the rest of treatments
  norm_treatment_means <- rowMeans(dplyr::select(metabolites_table, samples_norm))
  #print(norm_treatment_means)
  # Get the sds of the treatment used to normalize the rest of treatments
  #norm_treatment_stds <- apply(dplyr::select(metabolites_table, starts_with(treatment_to_normilize)), MARGIN = 1, FUN = sd)
  norm_treatment_stds <- apply(dplyr::select(metabolites_table, samples_norm), MARGIN = 1, FUN = sd)
  #print(norm_treatment_stds)
  
  for (metabolite in 1:nrow(metabolites_table_norm)) {
    #print(metabolite)
    for (cSample in 1:ncol(metabolites_table_norm)) {
      #print(metabolites_table_norm[metabolite, cSample]*1000)
      #metabolites_table_norm[metabolite, cSample] <- (metabolites_table_norm[metabolite, cSample]-norm_treatment_means[metabolite])/norm_treatment_stds[metabolite]
      metabolites_table_norm[metabolite, cSample] <- (metabolites_table_norm[metabolite, cSample]-norm_treatment_means[metabolite])
    }
  }
  
  #print(head(metabolites_table_norm))
  #rownames(metabolites_table_norm) <- metabolites_table$Metabolite
  metabolites_table_norm <- tibble::add_column(metabolites_table_norm, "Metabolite" = metabolites_table$Metabolite, .before = 1)
  return(metabolites_table_norm)
}

graph_metabolites <- function(metabolite_table_g, y1 = -10, y2= 10, dotsize = 0.5, binwidth = 0.5, ylab = "Log fold change"){
  dotplot <- ggplot(data = metabolite_table_g, aes(x = Metabolite, y = Change, fill = sample_type))
  dotplot + geom_dotplot(binaxis = 'y', dotsize = dotsize, stackdir='center', position=position_dodge(0.8), binwidth = binwidth) +
    labs(y= ylab, x = "Metabolites") +
    ylim(y1, y2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
}

transdata <- t(feature_table)
transdata <- as.data.frame(transdata)
colnames(transdata) <- transdata[2,]
transdata <- transdata[5:nrow(transdata), ]

#fiadata <- select(X20240822_FIA_Data_LC_Test_pos, c(2, 5:ncol(X20240822_FIA_Data_LC_Test_pos)))
#metNames <- fiadata$MetName
#fiadata <- fiadata[2:ncol(fiadata)]
#row.names(fiadata) <- metNames

new_sample_names <- paste(LC_metadata$syncom, LC_metadata$group, LC_metadata$replicate, sep = "_")

in_feat_table <- extract_features_comparison(feature_table = transdata, sig_features = c("L-Alanine[M+H]+", "L-Arginine[M+H]+", "L-Asparagine[M+H]+",
                                                                                          "L-Aspartic acid[M+H]+", "L-Cystine[M+H]+", "L-Glutamic acid[M+H]+",
                                                                                         "L-Glutamine[M+H]+", "L-Histidine[M+H]+", "L-Isoleucine[M+H]+",
                                                                                         "L-Leucine[M+H]+", "L-Lysine[M+H]+", "L-Methionine[M+H]+",
                                                                                         "L-Phenylalanine[M+H]+", "L-Proline[M+H]+", "L-Serine[M+H]+",
                                                                                         "L-Threonine[M+H]+", "L-Tryptophan[M+H]+", "L-Tyrosine[M+H]+",
                                                                                         "L-Valine[M+H]+"), columns_to_preserve = NULL)
in_feat_table <- in_feat_table

colnames(in_feat_table) <- c(new_sample_names, "Metabolite")


metab_flooding <- select(in_feat_table, contains("flooding"), "Metabolite")

metab_flooding_log <- log2_convert(metab_flooding)

norm_table_metab_flooding <- normalize_table_to_treatment(metabolites_table = metab_flooding_log, samples_norm = colnames(select(metab_flooding, contains("Control"))))

metab_flooding_g <- tidyr::gather(data = norm_table_metab_flooding, key = "Sample", value = "Change", colnames(norm_table_metab_flooding[,2:ncol(norm_table_metab_flooding)]))

metab_flooding_g$sample_type <- sapply(strsplit(metab_flooding_g$Sample, "_"), "[", 1)

graph_metabolites(metab_flooding_g, , y1 = -7.5, y2 = 3, dotsize = 25, binwidth = 0.01)



metab_scrapping <- select(in_feat_table, contains("scrapping"), "Metabolite")

metab_scrapping_log <- log2_convert(metab_scrapping)

norm_table_metab_scrapping <- normalize_table_to_treatment(metabolites_table = metab_scrapping_log, samples_norm = colnames(select(metab_scrapping, contains("Control"))))

metab_scrapping_g <- tidyr::gather(data = norm_table_metab_scrapping, key = "Sample", value = "Change", colnames(norm_table_metab_scrapping[,2:ncol(norm_table_metab_scrapping)]))

metab_scrapping_g$sample_type <- sapply(strsplit(metab_scrapping_g$Sample, "_"), "[", 1)

graph_metabolites(metab_scrapping_g, , y1 = -3, y2 = 7.5, dotsize = 25, binwidth = 0.01)

