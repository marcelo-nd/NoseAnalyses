# Load microbiome graph helper
source("https://raw.githubusercontent.com/marcelo-nd/microbiome-help/main/graphs.R")

otu_table_s1 <- read.csv("E:/1_NoseSynComProject/SequencingData/1_16S_First_Test_17_01_24/results/otu_table.csv", row.names=1)

otu_table_s2 <- read.csv("E:/1_NoseSynComProject/SequencingData/2_Test2_240124/results/otu_table.csv", row.names=1)

otu_table_s3 <- read.csv("E:/1_NoseSynComProject/SequencingData/3_SynComTest3/results/otu_table.csv", row.names=1)

otu_table_s4 <- read.csv("E:/1_NoseSynComProject/SequencingData/4_Test4-SynComNewPrimers/results/otu_table.csv", row.names=1)

#otu_table_s5 <- read.csv("E:/1_NoseSynComProject/SequencingData/Test4-SynComNewPrimers/results2/otu_table.csv", row.names=1)

otu_table_s5 <- read.csv("E:/1_NoseSynComProject/SequencingData/5_240305_Primer_27F-Test/results/otu_table.csv", row.names=1)

otu_table_s5 <- read.csv("E:/1_NoseSynComProject/SequencingData/5_240305_Primer_27F-Test/results/otu_table.csv", row.names=1)

otu_table_s7 <- read.csv("E:/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/results/otu_table.csv", row.names=1)

otu_table_s8 <- read.csv("E:/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/results_mirror/otu_table.csv", row.names=1)


# Raw results
barplot_from_feature_table(otu_table_s1)

barplot_from_feature_table(otu_table_s2)

barplot_from_feature_table(otu_table_s3)

barplot_from_feature_table(otu_table_s4)

barplot_from_feature_table(otu_table_s5)

barplot_from_feature_table(otu_table_s7)

barplot_from_feature_table(otu_table_s8)



sub_otutable <- function(otu_table, sample_indices, sample_names){
  # Select samples defined by user
  sub_otut <- dplyr::select(otu_table, any_of(sample_indices))
  # Remove rows that have only 0s
  sub_otut <- sub_otut[rowSums(sub_otut != 0) > 0, ]
  # if user provided sample names, apply them
  if (!is.null(sample_names) && length(sample_indices) == length(sample_names)) {
    colnames(sub_otut) <- sample_names
  }
  return(sub_otut)
}


### SynCom Inoculums

# 4 SynCom Inocs 27 F (from otu_table_s3)

inocs_27f <- sub_otutable(otu_table_s3, c(1:4), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_27f_plot <- barplot_from_feature_table(inocs_27f)

plot(inocs_27f_plot)

# 4 SynCom Inocs 27 F-II (from otu_table_s4)

inocs_27f2 <- sub_otutable(otu_table_s4, c(1:4), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_27f2_plot <- barplot_from_feature_table(inocs_27f2)

plot(inocs_27f2_plot)

# 4 SynCom Inocs V34 (Illumina Stlye) (from otu_table_s7)

inocs_v34 <- sub_otutable(otu_table_s7, c(5:8), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_v34_plot <- barplot_from_feature_table(inocs_v34)

plot(inocs_v34_plot)

# 4 SynCom Inocs rrna

inocs_rrna <- sub_otutable(otu_table_s8, c(5:8), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_rrna_plot <- barplot_from_feature_table(inocs_rrna)

plot(inocs_rrna_plot)


# Together
inocs_plot <- barplot_from_feature_tables(feature_tables = list(inocs_27f, inocs_27f2, inocs_v34, inocs_rrna), experiments_names = c("Inoc. 27F", "Inoc. 27FII", "Inoc. v34", "Inoc. rRNA"))

plot(inocs_plot)


### Mock DNA mixes

# Mock Mixes theoretical
otu_table_dna_mix <- readxl::read_excel("C:/Users/marce/OneDrive - UT Cloud/Mappe1.xlsx", sheet = "Tabelle2", range = "I16:M27")

species <- otu_table_dna_mix$species

otu_table_dna_mix <- as.data.frame(otu_table_dna_mix[,2:5])

rownames(otu_table_dna_mix) <- species

barplot_from_feature_table(otu_table_dna_mix)


# 4 Mock Mixes V34 (Illumina Stlye) (from otu_table_s7)

mocks_v34 <- sub_otutable(otu_table_s7, c(1:4), c("Mock 1", "Mock 2", "Mock 3", "Mock 4"))

mocks_v34_plot <- barplot_from_feature_table(mocks_v34)

plot(mocks_v34_plot)

# 4 Mock Mixes rrna primer

mocks_rrna <- sub_otutable(otu_table_s8, c(1:4), c("Mock 1", "Mock 2", "Mock 3", "Mock 4"))

mocks_rrna_plot <- barplot_from_feature_table(mocks_rrna)

plot(mocks_rrna_plot)

# Together
mocks_plot <- barplot_from_feature_tables(feature_tables = list(mocks_v34, mocks_rrna), experiments_names = c("Inoc. v34", "Inoc. rRNA"))

plot(mocks_plot)



### Zymo Community Standars


# Zymo Community Standars 27 F (from otu_table_s3)

zcs_27f <- sub_otutable(otu_table_s3, c(23:24), c("zcs1-27F", "zcs-27F_2"))

zcs_27f <- zcs_27f[1]

colnames(zcs_27f) <- "zcs"

zcs_27f_plot <- barplot_from_feature_table(zcs_27f)

plot(zcs_27f_plot)

# Zymo Community Standars 27 F-II (from otu_table_s4)

zcs_27f2 <- sub_otutable(otu_table_s4, c(13:14), c("zcs1-27FII", "zcs-27FII_2"))

zcs_27f2 <- zcs_27f2[1]

colnames(zcs_27f2) <- "zcs"

zcs_27f2_plot <- barplot_from_feature_table(zcs_27f2)

plot(zcs_27f2_plot)

# Zymo Community Standars V34 (Illumina Stlye) (from otu_table_s7)

zcs_v34 <- sub_otutable(otu_table_s7, c(9:10), c("zcs1-v34", "zcs-v34_2"))

zcs_v34 <- zcs_v34[1]

colnames(zcs_v34) <- "zcs"

zcs_v34_plot <- barplot_from_feature_table(zcs_v34)

plot(zcs_v34_plot)

# Zymo Community Standars rrna (from otu_table_s8)

zcs_rrna <- sub_otutable(otu_table_s8, c(9:10), c("zcs1-rrna", "zcs-rrna_2"))

zcs_rrna <- zcs_rrna[1]

colnames(zcs_rrna) <- "zcs"

zcs_rrna_plot <- barplot_from_feature_table(zcs_rrna)

plot(zcs_rrna_plot)

# Together
zcs_plot <- barplot_from_feature_tables(feature_tables = list(zcs_27f, zcs_27f2, zcs_v34, zcs_rrna), experiments_names = c("ZCS_27F", "ZCS_27FII", "ZCS_v34", "ZCS_rRNA"))

plot(zcs_plot)



### First 6 SynCom Experiment.

# 6 SynComs
syncomsp1 <- select(otu_table_s1, "barcode01", "barcode02", "barcode03")
barplot_from_feature_table(syncomsp1)
syncomsp2 <- select(otu_table_s2, "barcode01", "barcode02", "barcode03", "barcode04", "barcode05","barcode06","barcode07","barcode08","barcode09")
barplot_from_feature_table(syncomsp2)
syncomsp3 <- select(otu_table_s5, "barcode05", "barcode06", "barcode07", "barcode08", "barcode17", "barcode18")
barplot_from_feature_table(syncomsp3)


colnames(syncomsp1) <- c("SynCom1R1", "SynCom1R2", "SynCom1R3")
colnames(syncomsp2) <- c("SynCom2R1", "SynCom2R2", "SynCom2R3", "SynCom4R1", "SynCom4R2", "SynCom4R3", "SynCom6R1", "SynCom6R2", "SynCom6R3")
colnames(syncomsp3) <- c("SynCom3R1", "SynCom3R2", "SynCom3R3", "SynCom5R1", "SynCom5R2", "SynCom5R3")

syncomsp1$species <- rownames(syncomsp1)
syncomsp2$species <- rownames(syncomsp2)
syncomsp3$species <- rownames(syncomsp3)

syncomall <- dplyr::full_join(syncomsp1, syncomsp2, by="species")

syncomall <- dplyr::full_join(syncomall, syncomsp3, by="species")

row.names(syncomall) <- syncomall$species

syncomall <- select(syncomall, -c("species"))

syncomall2 <- syncomall %>% replace(is.na(.), 0)

syncomall2 <- syncomall2[rowSums(syncomall2 != 0) > 0, ]

syncomall3 <- syncomall2[!(row.names(syncomall2) %in% c("Escherichia coli")), ]

barplot_from_feature_table(syncomall3)
