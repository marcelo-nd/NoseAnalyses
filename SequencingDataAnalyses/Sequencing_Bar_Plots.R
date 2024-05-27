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

join_otu_tables <- function(otu_table, otu_table2){
  
}


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

plot(zcs_plot)

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



###############################################################

otu_table_s10 <- read.csv("F:/SequencingData/SynComTFBatch1and2_252524/results_long/otu_table.csv", row.names=1)

barplot_from_feature_table(otu_table_s10)

barplot_from_feature_table(otu_table_s10[1:61,])

otu_table_s10 <- otu_table_s10[1:61,]

# Select SynComs

sc1_ab <- sub_otutable(otu_table_s10, c(5:7), c("SC1_r1", "SC1_r2", "SC1_r3"))

barplot_from_feature_table(sc1_ab)

sc2_ab <- sub_otutable(otu_table_s10, c(4, 8, 9), c("SC2_r1", "SC2_r2", "SC2_r3"))

barplot_from_feature_table(sc2_ab)

sc22_ab <- sub_otutable(otu_table_s10, c(13:15), c("SC22_r1", "SC22_r2", "SC22_r3"))

barplot_from_feature_table(sc22_ab)

sc23_ab <- sub_otutable(otu_table_s10, c(16, 17, 18), c("SC23_r1", "SC23_r2", "SC23_r3"))

barplot_from_feature_table(sc23_ab)

sc24_ab <- sub_otutable(otu_table_s10, c(19, 20, 2), c("SC24_r1", "SC24_r2", "SC24_r3"))

barplot_from_feature_table(sc24_ab)

syncom_plot <- barplot_from_feature_tables(feature_tables = list(sc1_ab, sc2_ab, sc22_ab, sc23_ab, sc24_ab), experiments_names = c("SynCom1", "SynCom2", "SynCom3", "SynCom4", "SynCom5"))
plot(syncom_plot)


barplot_from_feature_tables <- function(feature_tables, experiments_names){
  # each experiment should have a separate table and should share samples and species
  for (table in seq(from = 1, to = length(feature_tables), by=1)) {
    #print(head(feature_tables))
    # copy feature table
    feature_table2 <- feature_tables[[table]]
    
    #print(head(feature_table2))
    
    # Remove columns (samples) with zero count
    if (ncol(feature_table2) > 1) {
      feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
    }
    
    # Generate a column with the names of ASVs/OTUs using rownames.
    feature_table2["species"] <- row.names(feature_table2)
    
    feature_table_g <- tidyr::gather(feature_table2, 1:(ncol(feature_table2) - 1) , key = "sample", value = "abundance")
    
    #print(experiments_names[table])
    
    feature_table_g$experiment <- experiments_names[table]
    
    if (table == 1) {
      exp_plot_table <- feature_table_g
    }else{
      exp_plot_table <- rbind(exp_plot_table, feature_table_g)
    }
  }
  
  #print(head(exp_plot_table))
  #print(head(feature_table2))
  
  # Keep order of samples in graph
  #exp_plot_table$sample <- factor(exp_plot_table$sample, levels = experiments_names)
  #print(head(exp_plot_table))
  # 
  color_palette <- get_palette()
  
  #exp_plot_table$experiment <- factor(exp_plot_table$experiment,levels=unique(as.character(exp_plot_table$experiment)))
  #exp_plot_table$species <- factor(exp_plot_table$species,levels = unique(as.character(exp_plot_table$species)))
  
  exp_plot_table <- exp_plot_table %>%
    dplyr::arrange(species)
  
  # Reorder factor
  #exp_plot_table$species <- forcats::fct_relevel(exp_plot_table$species, after = 0)
  #exp_plot_table$species <- forcats::fct_rev(exp_plot_table$species)
  
  print(exp_plot_table)
  
  otu_barplot <- ggplot(exp_plot_table) +
    geom_bar(aes(x = sample, y = abundance, fill = species),
             position = position_fill(),
             stat = "identity") + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::scale_fill_manual(values=color_palette) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   legend.title=ggplot2::element_text(size=14), 
                   legend.text=ggplot2::element_text(size=12)) +
    facet_grid(~experiment, scales = "free", space = "free") # this is to remove empty factors
  otu_barplot
  return(otu_barplot)
}

syncom_plot <- barplot_from_feature_tables(feature_tables = list(sc1_ab, sc2_ab, sc22_ab, sc23_ab, sc24_ab), experiments_names = c("SynCom1", "SynCom2", "SynCom22", "SynCom23", "SynCom24"))

plot(syncom_plot)
