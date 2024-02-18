table1 <- read.table("//wsl.localhost/Ubuntu/home/marcelo/results1/otu_table_sample1.tsv", sep=";", header = TRUE)

table2 <- read.table("//wsl.localhost/Ubuntu/home/marcelo/results2/otu_table_sample2.tsv", sep=";", header = TRUE)

table3 <- read.table("//wsl.localhost/Ubuntu/home/marcelo/results3/otu_table_sample3.tsv", sep=";", header = TRUE)

table4 <- read.table("//wsl.localhost/Ubuntu/home/marcelo/results4/otu_table_sample4.tsv", sep=";", header = TRUE)


main_otu_table <- table1

#table.list <- list(table2, table3, table4)
table.list <- list(table2)

for (table_number in seq(1, 1, 1)) {
  print(table_number)
  main_otu_table <- rbind(main_otu_table, table.list[[table_number]])
}

get_palette <- function(nColors = 50){
  return(c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2","brown1", "#CC79A7", "olivedrab3", "rosybrown",
           "darkorange3", "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
           "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
           "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
           "#e3ae78", "#a23f3f", "#290f76", "#ce7e00", "#386857",
           "#738564", "#e89d56", "#cd541d", "#1a3a46", "#ffe599",
           "#583E26", "#A78B71", "#F7C815", "#EC9704", "#9C4A1A",
           "firebrick2", "#C8D2D1", "#14471E", "#EE9B01", "#DA6A00",
           "#4B1E19", "#C0587E", "#FC8B5E", "#EA592A", "#FEF4C0")[1:nColors])
}

color_palette <- get_palette()

ggplot2::ggplot(main_otu_table, ggplot2::aes(x=sample, y=counts, fill=species)) + 
  ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
  ggplot2::scale_fill_manual(values=color_palette)
