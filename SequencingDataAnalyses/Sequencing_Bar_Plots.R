# Load microbiome graph helper
source("https://raw.githubusercontent.com/marcelo-nd/microbiome-help/main/graphs.R")

otu_table_s1 <- read.csv("E:/1_NoseSynComProject/SequencingData/1_16S_First_Test_17_01_24/results/otu_table.csv", row.names=1)

otu_table_s2 <- read.csv("E:/1_NoseSynComProject/SequencingData/2_Test2_240124/results/otu_table.csv", row.names=1)

otu_table_s3 <- read.csv("E:/1_NoseSynComProject/SequencingData/3_SynComTest3/results/otu_table.csv", row.names=1)

otu_table_s4 <- read.csv("E:/1_NoseSynComProject/SequencingData/Test4-SynComNewPrimers/results/otu_table.csv", row.names=1)

otu_table_s5 <- read.csv("E:/1_NoseSynComProject/SequencingData/Test4-SynComNewPrimers/results2/otu_table.csv", row.names=1)


barplot_from_feature_table(otu_table_s1)

barplot_from_feature_table(otu_table_s2)

# Only inocs
inoc <- otu_table_s2[,12:15]

inoc2 <- inoc[apply(inoc[,-1], 1, function(x) !all(x==0)),]

barplot_from_feature_table(inoc2)

#

barplot_from_feature_table(otu_table_s3)

barplot_from_feature_table(otu_table_s3[apply(otu_table_s3[,13:20], 1, function(x) !all(x==0)),][,13:20])

barplot_from_feature_table(otu_table_s3[apply(otu_table_s3[,5:12], 1, function(x) !all(x==0)),][,5:12])

barplot_from_feature_table(otu_table_s3[apply(otu_table_s3[,1:4], 1, function(x) !all(x==0)),][,1:4])

# Only inocs
row_sub1 = apply(otu_table_s3[,1:4], 1, function(row) all(row !=0 ))

barplot_from_feature_table(otu_table_s3[,1:4][row_sub1,])

#
barplot_from_feature_table(otu_table_s4)

row_sub2 = apply(otu_table_s4[,1:4], 1, function(row) all(row !=0 ))

barplot_from_feature_table(otu_table_s4[,1:4][row_sub2,])

#
barplot_from_feature_table(otu_table_s5[,1:4])


row_sub3 = apply(otu_table_s5[,1:4], 1, function(row) all(row !=0 ))

barplot_from_feature_table(otu_table_s5[,1:4][row_sub3,])






joined_tables <- dplyr::full_join(resultsotu_table, resultsotu_table_2)

joined_tables_e <- joined_tables[,-1]

rownames(joined_tables_e) <- joined_tables[,1]

joined_tables_e[is.na(joined_tables_e)] <- 0

barplot_from_feature_table(joined_tables_e)



standrd <- dplyr::select(joined_tables_e, "barcode04_1", "barcode10", "barcode11")

standrd <- standrd[rowSums(standrd != 0) > 0, ]

barplot_from_feature_table(standrd)



syncoms <- dplyr::select(joined_tables_e, "barcode01_2", "barcode02_2", "barcode03_1","barcode01", "barcode02", "barcode03", "barcode04", "barcode05", "barcode06", "barcode07", "barcode08", "barcode09")

syncoms <- syncoms[rowSums(syncoms != 0) > 0, ]

syncoms <- rbind(syncoms[1,], syncoms[4:7,])

barplot_from_feature_table(syncoms)


barplot_from_feature_table(otu_table)
