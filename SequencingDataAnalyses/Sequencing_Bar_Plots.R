# Load microbiome graph helper
source("https://raw.githubusercontent.com/marcelo-nd/microbiome-help/main/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")

otu_table_s1 <- read.csv("F:/SequencingData/Tests/1_16S_First_Test_170124/results/otu_table.csv", row.names=1)

otu_table_s2 <- read.csv("F:/SequencingData/Tests/2_Test2_240124/results/otu_table.csv", row.names=1)

otu_table_s3 <- read.csv("F:/SequencingData/Tests/3_SynComTest3_120224/results/otu_table.csv", row.names=1)

otu_table_s4 <- read.csv("F:/SequencingData/Tests/4_Test4-SynComNewPrimers_280224/results/otu_table.csv", row.names=1)

#otu_table_s4 <- read.csv("E:/1_NoseSynComProject/SequencingData/Test4-SynComNewPrimers/results2/otu_table.csv", row.names=1) # assigned with Emu DB

otu_table_s5 <- read.csv("F:/SequencingData/Tests/5_Primer_27F-Test_050324/results/otu_table.csv", row.names=1)

# 
otu_table_s6 <- read.csv("F:/SequencingData/Tests/6_Creat_rrna_070324/results_mirror/otu_table.csv", row.names=1)

# short reads of Run 7, processed with emu, LaCa 16S DB
otu_table_s7 <- read.csv("F:/SequencingData/Tests/7_PrimerConfTest_Mocks_110324/results/otu_table.csv", row.names=1)
# long reads of Run 7, processed with mirror
otu_table_s8 <- read.csv("F:/SequencingData/Tests/7_PrimerConfTest_Mocks_110324/results_mirror/otu_table.csv", row.names=1)

# long reads of Run 8, processed with mirror
otu_table_s9 <- read.csv("F:/SequencingData/Tests/8_confirmation_rrna_230324/results_mirror/otu_table.csv", row.names=1)

# Run 9, Mock DNA test using 16S barcoding Kit.
otu_table_s10 <- read.csv("F:/SequencingData/Tests/9_16STest_MockDNA_030424/results/otu_table.csv", row.names=1)



# Raw results
barplot_from_feature_table(otu_table_s1)

barplot_from_feature_table(otu_table_s2)

barplot_from_feature_table(otu_table_s3)

barplot_from_feature_table(otu_table_s4)

barplot_from_feature_table(otu_table_s5)
# short reads of Run 7, processed with emu
barplot_from_feature_table(otu_table_s7)
# long reads of Run 7, processed with mirror
barplot_from_feature_table(otu_table_s6)
# long reads of Run 7, processed with mirror
barplot_from_feature_table(otu_table_s8)
# long reads of Run 8, processed with mirror
barplot_from_feature_table(otu_table_s9)
# Run 9, Mock DNA test using 16S barcoding Kit.
barplot_from_feature_table(otu_table_s10)



### SynCom Inoculums
# Samples on Inoculums were sequenced only once with primers 27F, 27FII and V34. Twice with primer rrna.
# Since SynComs 1 and 2 are biological replicates, we really have duplicates for everything. and quatruplicates for rrna primers.

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
library(dplyr)
inocs_plot <- barplot_from_feature_tables(feature_tables = list(inocs_27f, inocs_27f2, inocs_v34, inocs_rrna), experiments_names = c("Inoc. 27F", "Inoc. 27FII", "Inoc. v34", "Inoc. rRNA"))

plot(inocs_plot)


### Mock DNA mixes
# Duplicates for V34 primers and triplicates for rrna.

# Mock Mixes theoretical
otu_table_dna_mix <- readxl::read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/DNA_Calculations.xlsx", sheet = "Nasal Mock Comms", range = "J18:N29")

species <- otu_table_dna_mix$species

otu_table_dna_mix <- as.data.frame(otu_table_dna_mix[,2:5])

rownames(otu_table_dna_mix) <- species

colnames(otu_table_dna_mix) <- c("Mock 1", "Mock 2", "Mock 3", "Mock 4")

barplot_from_feature_table(otu_table_dna_mix)


# 4 Mock Mixes 27F (Illumina Stlye) (from otu_table_s7)

mocks_27F <- sub_otutable(otu_table_s10, c(9, 12, 15, 2), c("Mock 1", "Mock 2", "Mock 3", "Mock 4"))

mocks_27F_plot <- barplot_from_feature_table(mocks_27F)

plot(mocks_27F_plot)


# 4 Mock Mixes V34 (Illumina Stlye) (from otu_table_s7)

mocks_v34 <- sub_otutable(otu_table_s7, c(1:4), c("Mock 1", "Mock 2", "Mock 3", "Mock 4"))

mocks_v34_plot <- barplot_from_feature_table(mocks_v34)

plot(mocks_v34_plot)


# 4 Mock Mixes rrna primer

mocks_rrna <- sub_otutable(otu_table_s8, c(1:4), c("Mock 1", "Mock 2", "Mock 3", "Mock 4"))

mocks_rrna_plot <- barplot_from_feature_table(mocks_rrna)

plot(mocks_rrna_plot)

# Together
mocks_plot <- barplot_from_feature_tables(feature_tables = list(otu_table_dna_mix, mocks_27F, mocks_v34, mocks_rrna), experiments_names = c("Theoretical", "Mock27F", "Mock v34", "Mock rRNA"))

mocks_plot <- barplot_from_feature_tables2(feature_tables = list(otu_table_dna_mix, mocks_27F, mocks_v34, mocks_rrna), experiments_names = c("Theoretical", "Mock27F", "Mock v34", "Mock rRNA"), shared_samples = TRUE)

plot(mocks_plot)


### Zymo Community Standards
# 27F: 10 seq replicates; 27FII: 2 replicates seq replicates; v34: 2 seq replicates; rrna: 6 seq replicates
# Theoretical
zcs_theo <- readxl::read_excel("C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/DNA_Calculations.xlsx", sheet = "zcs", range = "A2:B10")

species_zcs <- zcs_theo$species

zcs_theo <- as.data.frame(zcs_theo[,2])

rownames(zcs_theo) <- species_zcs

barplot_from_feature_table(zcs_theo)


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
zcs_plot <- barplot_from_feature_tables(feature_tables = list(zcs_theo, zcs_27f, zcs_27f2, zcs_v34, zcs_rrna), experiments_names = c("ZCS_theo", "ZCS_27F", "ZCS_27FII", "ZCS_v34", "ZCS_rRNA"))

zcs_plot2 <- barplot_from_feature_tables2(feature_tables = list(zcs_theo, zcs_27f, zcs_27f2, zcs_v34, zcs_rrna), experiments_names = c("ZCS_theo", "ZCS_27F", "ZCS_27FII", "ZCS_v34", "ZCS_rRNA"), shared_samples = TRUE)

plot(zcs_plot)

plot(zcs_plot2)

##### Mocks Triplicates
# Mock4
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

# Mock2
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

# Mock3
mock3_theo <- otu_table_dna_mix[3]
colnames(mock3_theo) <- "Mock3"
mock3_rrna1 <- otu_table_s8[3]
colnames(mock3_rrna1) <- "Mock3"
mock3_rrna2 <- otu_table_s9[3]
colnames(mock3_rrna2) <- "Mock3"
mock3_rrna3 <- otu_table_s9[8]
colnames(mock3_rrna3) <- "Mock3"

mock3_plot <- barplot_from_feature_tables(feature_tables = list(mock3_theo, mock3_rrna1, mock3_rrna2, mock3_rrna3), experiments_names = c("mock3_theo", "mock3_rRNA_1", "mock3_rRNA_2", "mock3_rRNA_3"))
plot(mock3_plot)

# Mock1
mock1_theo <- otu_table_dna_mix[1]
colnames(mock1_theo) <- "Mock1"
mock1_rrna1 <- otu_table_s8[1]
colnames(mock1_rrna1) <- "Mock1"
mock1_rrna2 <- otu_table_s9[1]
colnames(mock1_rrna2) <- "Mock1"
mock1_rrna3 <- otu_table_s9[6]
colnames(mock1_rrna3) <- "Mock1"

mock1_plot <- barplot_from_feature_tables(feature_tables = list(mock1_theo, mock1_rrna1, mock1_rrna2, mock1_rrna3), experiments_names = c("mock1_theo", "mock1_rRNA_1", "mock1_rRNA_2", "mock1_rRNA_3"))
plot(mock1_plot)

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




################

# long reads of Run 8, processed with mirror
otu_table_s9_2 <- read.csv("C:/Users/marce/Desktop/otu_table2.csv", row.names=1)

barplot_from_feature_table(otu_table_s9_2)


otu_table_s9_3 <- read.csv("C:/Users/marce/Desktop/otu_table3.csv", row.names=1)

barplot_from_feature_table(otu_table_s9_3)


###

otu_table_s8_emu <- read.csv("E:/1_NoseSynComProject/SequencingData/8_results_emu_RRNA_DB_no_repro/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s8_emu)



#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
source("C:/Users/marce/Documents/GitHub/microbiome-help/microbiomeGraphing.R")

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")


otu_table_s10 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/ReproData/SynComTFBatch1and2_250524/emu_results/otu_table.csv", row.names=1)

otu_table_s10 <- otu_table_s10[1:61,]

otu_table_s10 <- filter_otus_by_counts_col_counts(otu_table_s10, 10, 1)

barplot_from_feature_table(otu_table_s10)

# Select SynComs

sc1_ab <- sub_otutable(otu_table_s10, c(5:7), c("SC1_r1", "SC1_r2", "SC1_r3"))

barplot_from_feature_table(sc1_ab)

sc2_ab <- sub_otutable(otu_table_s10, c(4, 8, 9), c("SC2_r1", "SC2_r2", "SC2_r3"))

barplot_from_feature_table(sc2_ab)

sc22_ab <- sub_otutable(otu_table_s10, c(13:15), c("SC6_r1", "SC6_r2", "SC6_r3"))

barplot_from_feature_table(sc22_ab)

sc23_ab <- sub_otutable(otu_table_s10, c(16, 17, 18), c("SC7_r1", "SC7_r2", "SC7_r3"))

barplot_from_feature_table(sc23_ab)

sc24_ab <- sub_otutable(otu_table_s10, c(19, 20, 2), c("SC11_r1", "SC11_r2", "SC11_r3"))

barplot_from_feature_table(sc24_ab)

syncom_plot <- barplot_from_feature_tables(feature_tables = list(sc1_ab[1:3,], sc2_ab[1:3,], sc22_ab[1:3,], sc23_ab[1:6,], sc24_ab[1:4,]), experiments_names = c("SynCom1", "SynCom2", "SynCom06", "SynCom07", "SynCom11"), shared_samples =  FALSE)

plot(syncom_plot)

###############
otu_table_short <- read.csv("F:/SequencingData/NanoporeTech/rrnaTests/emu_results_short/otu_table.csv", row.names=1)
colnames(otu_table_short) <- c("Zymo Std R1","Zymo Std R2","Zymo Std R3", "NasalDNAMix R1", "NasalDNAMix R2",  "NasalDNAMix R3", "NasalCellMix R1", "NasalCellMix R2")
barplot_from_feature_table(otu_table_short)

otu_table_medium <- read.csv("F:/SequencingData/NanoporeTech/rrnaTests/emu_results_medium/otu_table.csv", row.names=1)
colnames(otu_table_medium) <- c("Zymo Std R1","Zymo Std R2","Zymo Std R3", "NasalDNAMix R1", "NasalDNAMix R2",  "NasalDNAMix R3", "NasalCellMix R1", "NasalCellMix R2")
barplot_from_feature_table(otu_table_medium)

otu_table_long <- read.csv("F:/SequencingData/NanoporeTech/rrnaTests/emu_results_long//otu_table.csv", row.names=1)
colnames(otu_table_long) <- c("Zymo Std R1","Zymo Std R2","Zymo Std R3", "NasalDNAMix R1", "NasalDNAMix R2",  "NasalDNAMix R3", "NasalCellMix R1", "NasalCellMix R2")
barplot_from_feature_table(otu_table_long)

otu_table_only_short <- read.csv("F:/SequencingData/NanoporeTech/rrnaTests/emu_results_only_short/otu_table.csv", row.names=1)
colnames(otu_table_only_short) <- c("Zymo Std R1","Zymo Std R2","Zymo Std R3", "NasalDNAMix R1", "NasalDNAMix R2",  "NasalDNAMix R3", "NasalCellMix R1", "NasalCellMix R2")
barplot_from_feature_table(otu_table_only_short)

barplot_from_feature_tables(list(otu_table_only_short, otu_table_short, otu_table_medium, otu_table_long), c("Short", ">= Short", ">= Medium", "<= Long"), shared_samples = TRUE)

barplot_from_feature_tables(list(otu_table_only_short[c(1, 4, 7)], otu_table_short[c(1, 4, 7)], otu_table_medium[c(1, 4, 7)], otu_table_long[c(1, 4, 7)]), c("Short", ">= Short", ">= Medium", "<= Long"), shared_samples = TRUE)


#####

otu_table_070624 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/Sequencing of Batches 1 and 2/SynComExp070624/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_070624)

otu_table_070624_filt <- filter_otus_by_counts_col_counts(otu_table_070624[1:19], min_count = 50,
                                                          col_number = 1)

barplot_from_feature_table(otu_table_070624_filt)



otu_table_070624_filt <- filter_otus_by_counts_col_counts(otu_table_070624, min_count = 50,
                                                          col_number = 1)

barplot_from_feature_table(otu_table_070624_filt[1:13,])



####### Karo 16S

otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/KaroSC_NasalSC/Karo/emu_results/otu_table.csv", row.names=1)

otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/20240829_Karo_PRimer_SoilSynCom2/no_sample/20240829_1253_MN45148_AVM806_08a69708/fastq_pass/emu_results/otu_table.csv", row.names=1)
# Raw results
barplot_from_feature_table(otu_table_s1[4:12])

####### Nose

otu_table_s1 <- read.csv("F:/SequencingData/SynCom100/TimePoints/emu_results/otu_table.csv", row.names=1)
# Raw results
barplot_from_feature_table(otu_table_s1[,1:5])


########
####### Karo 16S

otu_table_s1 <- read.csv("D:/SequencingData/Karo2/emu_results/otu_table.csv", row.names=1)
# Raw results
barplot_from_feature_table(dplyr::select(otu_table_s1, c(1,2,3,7,8,9)))

barplot_from_feature_table(dplyr::select(otu_table_s1, c(10:17)))

otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/NasalSC100_300924/no_sample_id/20240930_1731_MN45148_FBA32257_3c181886/fastq_pass/emu_results/otu_table.csv", row.names=1)
# Raw results
barplot_from_feature_table(otu_table_s1)

##### Sequencing 25.10.24
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SequencingRunBatch5/no_sample_id/20241025_1411_MN45148_FAZ32262_7e5d6e09/fastq_pass/nose/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s1, legend_pos = "bottom", legend_cols = 5)

gut_table <- otu_table_s1[, 15:23]

gut_table <- gut_table[apply(gut_table, 1, function(row) sum(row >= 10) >= 1), ]

barplot_from_feature_table(gut_table)

otu_table_s2 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SequencingRunBatch5/no_sample_id/20241025_1411_MN45148_FAZ32262_7e5d6e09/fastq_pass/nose/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s2)

##### Sequencing 28.10.24
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_SR2/no_sample_id/20241028_1159_MN45148_AVN094_845724a9/fastq_pass/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s1)

##### Sequencing 29.10.24
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SeqRun3_291024/no_sample_id/20241029_1750_MN45148_AUN642_9d368fc1/fastq_pass/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s1)

##### Sequencing 12.11.24
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SeqRun_SC100_121124/no_sample_id/20241112_1515_MN45148_auh943_d79d3d47/fastq_pass/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s1)

##### Sequencing 15.11.24
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_SeqRun_151124/no_sample_id/20241115_1500_MN45148_FBA32257_7de35303/fastq_pass/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

##### Sequencing 28.11.24
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_281124/no_sample_id/20241128_1507_MN45148_FAY35296_ff71f5b3/fastq_pass/emu_results/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

##### Karos Run
##### emuDB
otu_table_s1 <- read.csv("C:/Users/marce/Desktop/pass/emu_results_emuDB/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

##### karoDB
otu_table_s1 <- read.csv("C:/Users/marce/Desktop/pass/emu_results_karoDB/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

#### run 20.01.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_200125/no_sample_id/20250120_1209_MN45148_aui012_771ff743/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

#### run 23.01.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_230125/no_sample_id/20250123_1637_MN45148_AVM817_037deab6/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

#### run 29.01.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_290125/no_sample_id/20250129_1629_MN45148_AUI367_bf1e5858/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1[1:19], legend_pos = "right", legend_cols = 1)

# Duda SC12
otu_table_s1 <- read.csv("E:/SequencingData/SynCom100/duda/SC12/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# Duda SC13
otu_table_s1 <- read.csv("E:/SequencingData/SynCom100/duda/SC13/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 31.01.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_300125/no_sample_id/20250131_1506_MN45148_awo447_fd9b78e9/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 05.02.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_050224/no_sample_id/20250205_1052_MN45148_FBA33889_f4d732c7/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 07.02.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_070225/no_sample_id/20250207_1300_MN45148_FAY35296_64b2d4a6/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 13.02.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_130225/no_sample_id/20250213_1012_MN45148_AYU138_31d9206a/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 14.02.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_140225/no_sample_id/20250214_1618_MN45148_AYY707_b38a097e/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 18.02.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_180225/no_sample_id/20250218_1259_MN45148_AWS287_8c3aa7bd/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 21.02.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_200225/no_sample_id/20250220_1746_MN45148_FAZ32262_c1bd7e7d/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 07.03.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/NS_SC100_070325/no_sample_id/20250307_1222_MN45148_FBA33889_f5fe9dab/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1, legend_pos = "right", legend_cols = 1)

# run 07.04.25
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SCPlus_070425/no_sample_id/20250407_1236_MN45148_AYT875_52253316/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1[1:23], legend_pos = "right", legend_cols = 1)

#
otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100Plus_try2_160425/no_sample_id/20250416_1346_MN45148_ayz085_541b5a66/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1[1:6], legend_pos = "right", legend_cols = 1)

otu_table_s1 <- read.csv("D:/1_NoseSynComProject/SequencingData/OriginalRuns/SC100Plus_try2_160425/no_sample_id/20250416_1346_MN45148_ayz085_541b5a66/fastq_pass/emu_results/otu_table.csv", row.names=1)
barplot_from_feature_table(otu_table_s1[7:15], legend_pos = "right", legend_cols = 1)
