# Load microbiome graph helper
source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")

#source("https://raw.githubusercontent.com/marcelo-nd/microbiome-help/main/microbiomeGraphing.R")

#source("https://raw.githubusercontent.com/marcelo-nd/microbiome-help/main/otuTableWrangling.R")


# 27F
S27F <- read.csv("F:/SequencingData/NanoporeTech/27F/emu_results/otu_table.csv", row.names=1)
S27F <- S27F[1:15,]
S27F_filt <- filter_otus_by_counts_col_counts(S27F, min_count = 50, col_number = 3)

barplot_from_feature_table(S27F_filt)

# 27FII
S27FII <- read.csv("F:/SequencingData/NanoporeTech/27FII/emu_results/otu_table.csv", row.names=1)
S27FII <- S27FII[1:14,]
S27FII_filt <- filter_otus_by_counts_col_counts(S27FII, min_count = 20, col_number = 1)

barplot_from_feature_table(S27FII_filt)

# V34
V34 <- read.csv("F:/SequencingData/NanoporeTech/V34/emu_results/otu_table.csv", row.names=1)
V34 <- V34[1:18,]
V34_filt <- filter_otus_by_counts_col_counts(V34, min_count = 20, col_number = 1)

barplot_from_feature_table(V34_filt)

### rRNA
rrna <- read.csv("F:/SequencingData/NanoporeTech/rrna/emu_results/otu_table.csv", row.names=1)
rrna <- rrna[1:22,]
rrna_filt <- filter_otus_by_counts_col_counts(rrna, min_count = 50, col_number = 3)
barplot_from_feature_table(rrna_filt)

S27F$species <- rownames(S27F)
S27FII$species <- rownames(S27FII)
V34$species <- rownames(V34)
rrna$species <- rownames(rrna)


################### Zymo Community Standard ###################
zcs_theo <- readxl::read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/DNA_Calculations.xlsx", sheet = "zcs", range = "A2:B10")
species_zcs <- zcs_theo$species
zcs_theo <- as.data.frame(zcs_theo[,2])
rownames(zcs_theo) <- species_zcs
barplot_from_feature_table(zcs_theo)

zymo27F <- S27F_filt[c(1, 2, 3)]
colnames(zymo27F) <- c("Zymo_27F_R1", "Zymo_27F_R2", "Zymo_27F_R3")
zymo27FII <- S27FII_filt[c(1, 2)]
colnames(zymo27FII) <- c("Zymo_27FII_R1", "Zymo_27FII_R2")
zymoV34 <- V34_filt[c(1, 2)]
colnames(zymoV34) <- c("Zymo_V34_R1", "Zymo_V34_R2")
zymorrna <- rrna_filt[c(1, 2, 3)]
colnames(zymorrna) <- c("Zymo_rrna_R1", "Zymo_rrna_R2", "Zymo_rrna_R3")

barplot_from_feature_tables(feature_tables = list(zcs_theo, zymo27F, zymo27FII, zymoV34, zymorrna),
                            experiments_names = c("ZCS Theoretical", "ZCS 27F", "ZCS 27FII",
                                                  "ZCS V34", "ZCS rRNA"),
                            shared_samples = FALSE)


################### Nasal DNA Mix ###################
otu_table_dna_mix <- readxl::read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/DNA_Calculations.xlsx", sheet = "Nasal Mock Comms", range = "J18:N29")
species <- otu_table_dna_mix$species
otu_table_dna_mix <- as.data.frame(otu_table_dna_mix[,2:5])
rownames(otu_table_dna_mix) <- species
colnames(otu_table_dna_mix) <- c("Mock 1", "Mock 2", "Mock 3", "Mock 4")
barplot_from_feature_table(otu_table_dna_mix)
NasalDNATheo <- otu_table_dna_mix[4]
colnames(NasalDNATheo) <- c("NasalDNA_Theoretical")

nasaldnaMock27F <- S27F_filt[c(4, 5, 6)]
colnames(nasaldnaMock27F) <- c("NasalDNA_27F_R1", "NasalDNA_27F_R2", "NasalDNA_27F_R3")
nasaldnaMock27FII <- S27FII_filt[c(3, 4, 5)]
colnames(nasaldnaMock27FII) <- c("NasalDNA_27FII_R1", "NasalDNA_27FII_R2", "NasalDNA_27FII_R3")
nasaldnaMockV34 <- V34_filt[c(3, 4, 5)]
colnames(nasaldnaMockV34) <- c("NasalDNA_V34_R1", "NasalDNA_V34_R2", "NasalDNA_V34_R3")
nasaldnaMockrrna <- rrna_filt[c(4, 5, 6)]
colnames(nasaldnaMockrrna) <- c("NasalDNA_rrna_R1", "NasalDNA_rrna_R2", "NasalDNA_rrna_R3")

barplot_from_feature_tables(feature_tables = list(NasalDNATheo, nasaldnaMock27F,
                                                  nasaldnaMock27FII, nasaldnaMockV34,
                                                  nasaldnaMockrrna),
                            experiments_names = c("Nasal DNA Theo.", "Nasal DNA 27F", "Nasal DNA 27FII",
                                                  "Nasal DNA V34", "Nasal DNA rRNA"),
                            shared_samples = FALSE)


################### Nasal Cell Mix ###################
nasalcellMock27F <- S27F_filt[c(7, 8, 9)]
colnames(nasalcellMock27F) <- c("NasalCELL_27F_R1", "NasalCELL_27F_R2", "NasalCELL_27F_R3")
nasalcellMock27FII <- S27FII_filt[c(6, 7, 8)]
colnames(nasalcellMock27FII) <- c("NasalCELL_27FII_R1", "NasalCELL_27FII_R2", "NasalCELL_27FII_R3")
nasalcellMockV34 <- V34_filt[c(6, 7, 8)]
colnames(nasalcellMockV34) <- c("NasalCELL_V34_R1", "NasalCELL_V34_R2", "NasalCELL_V34_R3")
nasalcellMockrrna <- rrna_filt[c(7, 8, 9)]
colnames(nasalcellMockrrna) <- c("NasalCELL_rrna_R1", "NasalCELL_rrna_R2", "NasalCELL_rrna_R3")

barplot_from_feature_tables(feature_tables = list(nasalcellMock27F,
                                                  nasalcellMock27FII, nasalcellMockV34,
                                                  nasalcellMockrrna),
                            experiments_names = c("Nasal CELL 27F", "Nasal CELL 27FII",
                                                  "Nasal CELL V34", "Nasal CELL rRNA"),
                            shared_samples = FALSE)
