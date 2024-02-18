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
