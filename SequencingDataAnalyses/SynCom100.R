# Load helping functions
source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")

# Screening

otu_table_screening <- read.csv("F:/SequencingData/SynCom100/Screening/emu_results/otu_table.csv", row.names=1)

colnames(otu_table_screening) <- c("SC1","SC2","SC3", "SC4", "SC5", "SC6", "SC7", "SC8",
                                   "SC9", "SC10", "SC11", "SC12", "SC13", "SC14", "SC15",
                                   "SC16", "SC17", "SC18", "SC19", "SC20", "SC21", "SC22",
                                   "SC23", "SC24")

otu_table_screening_filt <- filter_otus_by_counts_col_counts(otu_table_screening, min_count = 50,
                                                          col_number = 1)

barplot_from_feature_table(otu_table_screening_filt)
