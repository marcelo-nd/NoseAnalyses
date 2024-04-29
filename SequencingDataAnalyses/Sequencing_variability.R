################ Plots comparing technical variability

# Join two tables

# Load microbiome graph helper
source("https://raw.githubusercontent.com/marcelo-nd/microbiome-help/main/graphs.R")

otu_table_s1 <- read.csv("E:/1_NoseSynComProject/SequencingData/1_16S_First_Test_170124/results/otu_table.csv", row.names=1)

otu_table_s2 <- read.csv("E:/1_NoseSynComProject/SequencingData/2_Test2_240124/results/otu_table.csv", row.names=1)

otu_table_s3 <- read.csv("E:/1_NoseSynComProject/SequencingData/3_SynComTest3_120224/results/otu_table.csv", row.names=1)

otu_table_s4 <- read.csv("E:/1_NoseSynComProject/SequencingData/4_Test4-SynComNewPrimers_280224/results/otu_table.csv", row.names=1)

#otu_table_s4 <- read.csv("E:/1_NoseSynComProject/SequencingData/Test4-SynComNewPrimers/results2/otu_table.csv", row.names=1) # assigned with Emu DB

otu_table_s5 <- read.csv("E:/1_NoseSynComProject/SequencingData/5_Primer_27F-Test_050324/results/otu_table.csv", row.names=1)

# 
otu_table_s6 <- read.csv("E:/1_NoseSynComProject/SequencingData/6_Creat_rrna_070324/results_mirror/otu_table.csv", row.names=1)

# short reads of Run 7, processed with emu, LaCa 16S DB
otu_table_s7 <- read.csv("E:/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks_110324/results/otu_table.csv", row.names=1)
# long reads of Run 7, processed with mirror
otu_table_s8 <- read.csv("E:/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks_110324/results_mirror/otu_table.csv", row.names=1)

# long reads of Run 8, processed with mirror
otu_table_s9 <- read.csv("E:/1_NoseSynComProject/SequencingData/8_confirmation_rrna_230324/results_mirror/otu_table.csv", row.names=1)

# Run 9, Mock DNA test using 16S barcoding Kit.
otu_table_s10 <- read.csv("E:/1_NoseSynComProject/SequencingData/9_16STest_MockDNA_030424/results/otu_table.csv", row.names=1)


################################### Create datraframe with replicates and labeling for PCA

######################################################################## Zymo Community Standard

# Theoretical
zcs_theo <- readxl::read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/DNA_Calculations.xlsx", sheet = "zcs", range = "A2:B10")

species_zcs <- zcs_theo$species

zcs_theo <- as.data.frame(zcs_theo[,2])

rownames(zcs_theo) <- species_zcs

barplot_from_feature_table(zcs_theo)


# ZCS 27f
zcs_27f1 <- otu_table_s10[24]
colnames(zcs_27f1) <- "zcs_27f1"
zcs_27f1$species <- row.names(zcs_27f1)

zcs_27f2 <- otu_table_s5[19]
colnames(zcs_27f2) <- "zcs_27f2"
zcs_27f2$species <- row.names(zcs_27f2)

zcs_27f3 <- otu_table_s2[11]
colnames(zcs_27f3) <- "zcs_27f3"
zcs_27f3$species <- row.names(zcs_27f3)

zcs_27f_pca <- full_join(zcs_27f1, zcs_27f2, by = "species")
zcs_27f_pca <- full_join(zcs_27f_pca, zcs_27f3, by = "species")

row.names(zcs_27f_pca) <- zcs_27f_pca$species

zcs_27f_pca <- zcs_27f_pca[, -2]

zcs_27f_pca[is.na(zcs_27f_pca)] <- 0

barplot_from_feature_table(zcs_27f_pca)


# ZCS 27f2
zcs_27f2_1 <- otu_table_s4[13]
colnames(zcs_27f2_1) <- "zcs_27f2_1"
zcs_27f2_1$species <- row.names(zcs_27f2_1)

zcs_27f2_2 <- otu_table_s4[14]
colnames(zcs_27f2_2) <- "zcs_27f2_2"
zcs_27f2_2$species <- row.names(zcs_27f2_2)

zcs_27f2_pca <- full_join(zcs_27f2_1, zcs_27f2_2, by = "species")

row.names(zcs_27f2_pca) <- zcs_27f2_pca$species

zcs_27f2_pca <- zcs_27f2_pca[, -2]

zcs_27f2_pca[is.na(zcs_27f2_pca)] <- 0

barplot_from_feature_table(zcs_27f2_pca)


# ZCS v34
zcs_v34_1 <- otu_table_s7[9]
colnames(zcs_v34_1) <- "zcs_v34_1"
zcs_v34_1$species <- row.names(zcs_v34_1)

zcs_v34_2 <- otu_table_s7[10]
colnames(zcs_v34_2) <- "zcs_v34_2"
zcs_v34_2$species <- row.names(zcs_v34_2)

zcs_v34_pca <- full_join(zcs_v34_1, zcs_v34_2, by = "species")

row.names(zcs_v34_pca) <- zcs_v34_pca$species

zcs_v34_pca <- zcs_v34_pca[, -2]

zcs_v34_pca[is.na(zcs_v34_pca)] <- 0

barplot_from_feature_table(zcs_v34_pca)


# 27 rrna

##### ZCS rrna Triplicates
zcs_rrna1 <- otu_table_s8[9]
colnames(zcs_rrna1) <- "zcs_rrna_1"
#zcs_rrna2 <- otu_table_s8[10]
zcs_rrna2 <- otu_table_s6[1]
colnames(zcs_rrna2) <- "zcs_rrna_2"
zcs_rrna3 <- otu_table_s9[10]
colnames(zcs_rrna3) <- "zcs_rrna_3"

zcs_rrna1$species <- row.names(zcs_rrna1)
zcs_rrna2$species <- row.names(zcs_rrna2)
zcs_rrna3$species <- row.names(zcs_rrna3)

zcs_rrna_pca <- full_join(zcs_rrna1, zcs_rrna2, by = "species")

zcs_rrna_pca <- full_join(zcs_rrna_pca, zcs_rrna3, by = "species")

row.names(zcs_rrna_pca) <- zcs_rrna_pca$species

zcs_rrna_pca <- zcs_rrna_pca[, -2]

zcs_rrna_pca[is.na(zcs_rrna_pca)] <- 0

barplot_from_feature_table(zcs_rrna_pca)


######################################################################## Nasal Bacterial Mixes

# 27F
inocs_27f <- sub_otutable(otu_table_s3, c(1:4), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_27f_pca <- inocs_27f[, c(1,3)]

colnames(inocs_27f_pca) <- c("Inoc_27f_1", "Inoc_27f_2")

inocs_27f_pca$species <- row.names(inocs_27f_pca)

barplot_from_feature_table(inocs_27f_pca)

# 
inocs_27f2 <- sub_otutable(otu_table_s4, c(1:4), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_27f2_pca <- inocs_27f2[, c(1,3)]

colnames(inocs_27f2_pca) <- c("Inoc_27f2_1", "Inoc_27f2_2")

inocs_27f2_pca$species <- row.names(inocs_27f2_pca)

barplot_from_feature_table(inocs_27f2_pca)

#
inocs_v34 <- sub_otutable(otu_table_s7, c(5:8), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_v34_pca <- inocs_v34[, c(1,3)]

colnames(inocs_v34_pca) <- c("Inoc_v34_1", "Inoc_v34_2")

inocs_v34_pca$species <- row.names(inocs_v34_pca)

barplot_from_feature_table(inocs_v34_pca)

#
inocs_rrna <- sub_otutable(otu_table_s8, c(5:8), c("SynCom1 T0 0,1", "SynCom1 T0 0,01", "SynCom2 T0 0,1", "SynCom2 T0 0,01"))

inocs_rrna_pca <- inocs_rrna[, c(1,3)]
colnames(inocs_rrna_pca) <- c("Inoc_rrna_1", "Inoc_rrna_2")
inocs_rrna_pca$species <- row.names(inocs_rrna_pca)

inocs_rrna2 <- otu_table_s6[2]
colnames(inocs_rrna2) <- c("Inoc_rrna_3")
inocs_rrna2$species <- row.names(inocs_rrna2)

inocs_rrna_pca <- full_join(inocs_rrna_pca, inocs_rrna2, by = "species")

#inocs_rrna_pca <- inocs_rrna_pca[, -3]

#inocs_rrna_pca[is.na(inocs_rrna_pca)] <- 0

#barplot_from_feature_table(inocs_rrna_pca)


#barplot_from_feature_tables(feature_tables = list(inocs_27f_pca, inocs_27f2_pca, inocs_v34_pca, inocs_rrna_pca), experiments_names = c("Inoc. 27F", "Inoc. 27FII", "Inoc. v34", "Inoc. rRNA"))

inocs_df <- full_join(inocs_27f_pca, inocs_27f2_pca, by = "species")

inocs_df <- full_join(inocs_df, inocs_v34_pca, by = "species")

inocs_df <- full_join(inocs_df, inocs_rrna_pca, by = "species")

row.names(inocs_df) <- inocs_df$species

inocs_df <- inocs_df[, -3]

inocs_df[is.na(inocs_df)] <- 0

barplot_from_feature_table(inocs_df)

######################################################################## Mock DNA mixes
# Mock Mixes theoretical
otu_table_dna_mix <- readxl::read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/DNA_Calculations.xlsx", sheet = "Nasal Mock Comms", range = "J18:N29")
species <- otu_table_dna_mix$species

otu_table_dna_mix <- as.data.frame(otu_table_dna_mix[,2:5])

rownames(otu_table_dna_mix) <- species

colnames(otu_table_dna_mix) <- c("Mock 1", "Mock 2", "Mock 3", "Mock 4")

# Mock4 v34
mock4_theo <- otu_table_dna_mix[4]
colnames(mock4_theo) <- "Mock4"
mock4_v34_1 <- otu_table_s7[1]
colnames(mock4_v34_1) <- "Mock4"

mock4_v34_plot <- barplot_from_feature_tables(feature_tables = list(mock4_theo, mock4_v34_1), experiments_names = c("mock4_theo", "mock4_v34_1"))
plot(mock4_v34_plot)

# Mock2 V34
mock2_theo <- otu_table_dna_mix[2]
colnames(mock2_theo) <- "Mock2"
mock2_v34_1 <- otu_table_s8[2]
colnames(mock2_v34_1) <- "Mock2"

mock2_v34_plot <- barplot_from_feature_tables(feature_tables = list(mock2_theo, mock2_v34_1), experiments_names = c("mock2_theo", "mock2_v34_1"))
plot(mock2_v34_plot)


# Mock4 v34
mock4_theo <- otu_table_dna_mix[4]
colnames(mock4_theo) <- "Mock4"
mock4_v34_1 <- otu_table_s7[1]
colnames(mock4_v34_1) <- "Mock4"

mock4_v34_plot <- barplot_from_feature_tables(feature_tables = list(mock4_theo, mock4_v34_1), experiments_names = c("mock4_theo", "mock4_v34_1"))
plot(mock4_v34_plot)


# Mock4 rrna
mock4_theo <- otu_table_dna_mix[4]
colnames(mock4_theo) <- "Mock4"
mock4_rrna1 <- otu_table_s8[4]
colnames(mock4_rrna1) <- "Mock4"
mock4_rrna2 <- otu_table_s9[4]
colnames(mock4_rrna2) <- "Mock4"
mock4_rrna3 <- otu_table_s9[9]
colnames(mock4_rrna3) <- "Mock4"

mock4_plot <- barplot_from_feature_tables(feature_tables = list(mock4_theo, mock4_rrna1, mock4_rrna2, mock4_rrna3), experiments_names = c("mock4_theo", "mock4_rRNA_1", "mock4_rRNA_2", "mock4_rRNA_3"))
plot(mock4_plot)

# Mock2 rrna
mock2_theo <- otu_table_dna_mix[2]
colnames(mock2_theo) <- "Mock2"
mock2_rrna1 <- otu_table_s8[2]
colnames(mock2_rrna1) <- "Mock2"
mock2_rrna2 <- otu_table_s9[2]
colnames(mock2_rrna2) <- "Mock2"
mock2_rrna3 <- otu_table_s9[7]
colnames(mock2_rrna3) <- "Mock2"

mock2_plot <- barplot_from_feature_tables(feature_tables = list(mock2_theo, mock2_rrna1, mock2_rrna2, mock2_rrna3), experiments_names = c("mock2_theo", "mock2_rRNA_1", "mock2_rRNA_2", "mock2_rRNA_3"))
plot(mock2_plot)


# Mock4 27F
mock4_theo <- otu_table_dna_mix[4]
colnames(mock4_theo) <- "Mock4"
mock4_rrna1 <- otu_table_s8[4]
colnames(mock4_rrna1) <- "Mock4"
mock4_rrna2 <- otu_table_s9[4]
colnames(mock4_rrna2) <- "Mock4"
mock4_rrna3 <- otu_table_s9[9]
colnames(mock4_rrna3) <- "Mock4"

mock4_plot <- barplot_from_feature_tables(feature_tables = list(mock4_theo, mock4_rrna1, mock4_rrna2, mock4_rrna3), experiments_names = c("mock4_theo", "mock4_rRNA_1", "mock4_rRNA_2", "mock4_rRNA_3"))
plot(mock4_plot)

# Mock2 27F
mock2_theo <- otu_table_dna_mix[2]
colnames(mock2_theo) <- "Mock2"
mock2_rrna1 <- otu_table_s8[2]
colnames(mock2_rrna1) <- "Mock2"
mock2_rrna2 <- otu_table_s9[2]
colnames(mock2_rrna2) <- "Mock2"
mock2_rrna3 <- otu_table_s9[7]
colnames(mock2_rrna3) <- "Mock2"

mock2_plot <- barplot_from_feature_tables(feature_tables = list(mock2_theo, mock2_rrna1, mock2_rrna2, mock2_rrna3), experiments_names = c("mock2_theo", "mock2_rRNA_1", "mock2_rRNA_2", "mock2_rRNA_3"))
plot(mock2_plot)
