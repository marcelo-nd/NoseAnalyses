library(readxl)

library(ggplot2)

library(dplyr)

library(tidyverse)

library(scales) 

# Read Data

wd <- "D:/3/Still important/8_T-Metabolomics/GramPositiveTest_10_07_23/"

ratioList <- read_excel(paste0(wd,  "20230718_ratioList_g_post_test_ratios.xlsx"), sheet = "Ratios")

metabolites_ratios <- dplyr::filter(ratioList, Add == 1)

write.table(metabolites_ratios, file = "C:/Users/Marcelo/OneDrive - UT Cloud/1_Postdoc Tü/Sci/NoseSynComProject/8_T-Metabolomics/metabolites_ratios.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";")

##### Helper Functions
### 
normalize_table_to_treatment <- function(metabolites_table, treatment_to_normilize){
  #metabolites_table_norm <- metabolites_table %>% dplyr::select(-Metabolite) %>% t()
  metabolites_table_norm <- metabolites_table %>% dplyr::select(-Metabolite)
  
  # Get the means of the treatment used to normalize the rest of treatments
  norm_treatment_means <- rowMeans(dplyr::select(metabolites_table, starts_with(treatment_to_normilize)))
  #print(norm_treatment_means)
  # Get the sds of the treatment used to normalize the rest of treatments
  #norm_treatment_stds <- apply(dplyr::select(metabolites_table, starts_with(treatment_to_normilize)), MARGIN = 1, FUN = sd)
  norm_treatment_stds <- apply(dplyr::select(metabolites_table, starts_with(treatment_to_normilize)), MARGIN = 1, FUN = sd)
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

###
log2_table <- function(metabolites_table){
  metabolites_table_log2 <- metabolites_table %>% dplyr::select(-Metabolite) %>% log2()
  #metabolites_table_log2 <- metabolites_table %>% log2()
  #rownames(metabolites_table_log2) <- metabolites_table$Metabolite
  metabolites_table_log2 <- tibble::add_column(metabolites_table_log2, "Metabolite" = metabolites_table$Metabolite, .before = 1)
  return(metabolites_table_log2)
}

# Gathered samples
gather_samples <- function(metabolites_table, sample.names = c(), sample_type = NULL){
  #Metabolite <- rownames(metabolites_table)
  #metabolites_table_g <- cbind(Metabolite, metabolites_table)
  #print(head(metabolites_table_g))
  if (length(sample.names > 0)) {
    metabolites_table_g <- dplyr::select(metabolites_table, c("Metabolite", sample.names))
  }
  #print(head(metabolites_table_g))
  
  metabolites_table_g <- tidyr::gather(data = metabolites_table_g, key = "Sample", value = "Change", sample.names)
  if (!is.null(sample_type)) {
    metabolites_table_g$sample_type <- sample_type
  }
  
  return(metabolites_table_g)
}

###
graph_metabolites <- function(metabolite_table_g, y1 = -10, y2= 10, dotsize = 0.5, binwidth = 0.5, ylab = "12C/13C ratios fold change (log2)"){
  dotplot <- ggplot(data = metabolite_table_g, aes(x = Metabolite, y = Change, fill = sample_type))
  dotplot + geom_dotplot(binaxis = 'y', dotsize = dotsize, stackdir='center', position=position_dodge(0.8), binwidth = binwidth) +
    labs(y= ylab, x = "Metabolites") +
    ylim(y1, y2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
}

# E. coli
eco_metabolites <- dplyr::select(metabolites_ratios, Metabolite, Ecobb1, Ecobb2, Ecobb3, Ecoq1, Ecoq2, Ecoq3)

eco_metabolites_log2 <- log2_table(eco_metabolites)

eco_metabolites_std <- normalize_table_to_treatment(eco_metabolites_log2, "Ecoq")

Ecoq <- gather_samples(eco_metabolites_std, sample.names = c("Ecoq1", "Ecoq2", "Ecoq3"), sample_type = "E. coli qs")

Ecobb <- gather_samples(eco_metabolites_std, sample.names = c("Ecobb1", "Ecobb2", "Ecobb3"), sample_type = "E. coli qs + bb")

Eco_data <- rbind(Ecobb, Ecoq)

graph_metabolites(Eco_data, y1 = -5, y2 = 5, dotsize = 10, binwidth = 0.01)


# S. aureus
sau_metabolites <- dplyr::select(metabolites_ratios, Metabolite, Saubb1, Saubb2, Saubb3, Sauq1, Sauq2, Sauq3)

sau_metabolites_log2 <- log2_table(sau_metabolites)

sau_metabolites_std <- normalize_table_to_treatment(sau_metabolites_log2, "Sauq")

Sauq <- gather_samples(sau_metabolites_std, c("Sauq1", "Sauq2", "Sauq3"), sample_type = "S. aureus qs")

Saubb <- gather_samples(sau_metabolites_std, c("Saubb1", "Saubb2", "Saubb3"), sample_type = "S. aureus qs + bb")

Sau_data <- rbind(Saubb, Sauq)

graph_metabolites(Sau_data, y1 = -6, y2 = 5, dotsize = 10, binwidth = 0.01)

####


# Calculate Relative Standard deviation (RSD). This is calculated only with 12C signals.

x12CList <- read_excel("D:/3/Still important/8_T-Metabolomics/GramPositiveTest_10_07_23/20230718_12C_g_post_test_12c.xlsx",
                       sheet = "12C")

metabolites12c <- dplyr::filter(x12CList, Add == 1)

metabolites12c <- dplyr::select(metabolites12c, Metabolite, Ecobb1, Ecobb2, Ecobb3, Ecoq1, Ecoq2, Ecoq3, Saubb1, Saubb2, Saubb3, Sauq1, Sauq2, Sauq3)

metabolites12c <- metabolites12c[order(metabolites12c$Metabolite),]


Sauq_12c <- select(metabolites12c, c("Sauq1", "Sauq2", "Sauq3"))

Sauq_means <- apply(Sauq_12c, FUN = mean, MARGIN = 1)

Sauq_stds <- apply(Sauq_12c, FUN = sd, MARGIN = 1)

Sauq_rsds <- (Sauq_stds/Sauq_means)*100


Saubb_12c <- select(metabolites12c, c("Saubb1", "Saubb2", "Saubb3"))

Saubb_means <- apply(Saubb_12c, FUN = mean, MARGIN = 1)

Saubb_stds <- apply(Saubb_12c, FUN = sd, MARGIN = 1)

Saubb_rsds <- (Saubb_stds/Saubb_means)*100


Ecoq_12c <- select(metabolites12c, sample.names = c("Ecoq1", "Ecoq2", "Ecoq3"))

Ecoq_means <- apply(Ecoq_12c, FUN = mean, MARGIN = 1)

Ecoq_stds <- apply(Ecoq_12c, FUN = sd, MARGIN = 1)

Ecoq_rsds <- (Ecoq_stds/Ecoq_means)*100


Ecobb_12c <- select(metabolites12c, sample.names = c("Ecobb1", "Ecobb2", "Ecobb3"))

Ecobb_means <- apply(Ecobb_12c, FUN = mean, MARGIN = 1)

Ecobb_stds <- apply(Ecobb_12c, FUN = sd, MARGIN = 1)

Ecobb_rsds <- (Ecobb_stds/Ecobb_means)*100


graph_metabolites <- function(metabolite_table_g, y1 = -10, y2= 10, dotsize = 0.5, binwidth = 0.5){
  dotplot <- ggplot(data = metabolite_table_g, aes(x = Metabolite, y = Change, fill = sample_type))
  dotplot2 <- dotplot + geom_dotplot(binaxis = 'y', dotsize = dotsize, stackdir='center', position=position_dodge(0.8), binwidth = binwidth) +
    labs(y= "12C/13C ratios fold change (log2)", x = "Metabolites", size = 12) +
    ylim(y1, y2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), axis.text=element_text(size=11))
  
  return(dotplot2)
}



##### S. aureus graph

sau_graph <- graph_metabolites(Sau_data, y1 = -6, y2 = 5, dotsize = 10, binwidth = 0.01)

Sauq_rsds <- paste0(as.character(round(Sauq_rsds, digits = 2)), "%")

Saubb_rsds <- paste0(as.character(round(Saubb_rsds, digits = 2)), "%")

sau_graph +
  annotate("text",
           x=1:19, y=3.5, label= Sauq_rsds, size = 3, colour = "red") +
  annotate("text",
           x=1:19, y=3, label= Saubb_rsds, size = 3, colour = "blue")


##### E. coli graph

eco_graph <- graph_metabolites(Eco_data, y1 = -6, y2 = 5, dotsize = 10, binwidth = 0.01)

Ecoq_rsds <- paste0(as.character(round(Ecoq_rsds, digits = 2)), "%")

Ecobb_rsds <- paste0(as.character(round(Ecobb_rsds, digits = 2)), "%")

eco_graph +
  annotate("text",
           x=1:19, y=3.5, label= Ecoq_rsds, size = 3, colour = "red") +
  annotate("text",
           x=1:19, y=3, label= Ecobb_rsds, size = 3, colour = "blue")








############################


sau_rdsd <- c()
for (value in 1:19) {
  sau_rdsd <- c(sau_rdsd, paste(as.character(round(Sauq_rsds[value], digits = 2)), "%", "\n", as.character(round(Saubb_rsds[value], digits = 2)), "%"))
}

print(sau_rdsd)

#sep="\n",

eco_rdsd <- c()
for (value in 1:19) {
  eco_rdsd <- c(eco_rdsd, paste(as.character(round(Ecoq_rsds[value], digits = 2)), "%", "\n", as.character(round(Ecobb_rsds[value], digits = 2)), "%"))
}

print(eco_rdsd)

