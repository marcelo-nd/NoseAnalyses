# Bioconductor Packages

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install()
}


if (!"Biostrings" %in% installed.packages()) {
  BiocManager::install(c("Biostrings"))
}
library("Biostrings")

# Read fasta file
data_fasta <- readDNAStringSet("F:/16Sdatabases/RRN_db/species_taxid.fasta")

# Read taxonomy table
taxonomy <- read_delim("F:/16Sdatabases/RRN_db/taxonomy.tsv",
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)

extract_seqs_from_fasta_name <- function(fasta_file, species_names){
  fasta_results <- DNAStringSet() # Empty fasta file to store results
  
  # Iterate over original fasta file to extract from
  for(fasta_seq in 1:length(data_fasta)){
    # Iterate over species names
    for (specie in species_names) {
      current_seq = data_fasta[fasta_seq]
      # Look for name pattern string in all sequences names
      if(grepl(pattern = specie, x = current_seq@ranges@NAMES, fixed = TRUE)){
        print(paste("Found sequence for species: ", specie))
        fasta_results = append(x=fasta_results, current_seq) # append extracted sequences
      }
    }
  }
  print("Sequences extraction complete")
  
  # Check if all species are in extracted fasta file.
  found_species_fasta <- c()
  # Iterate over results
  for (fasta_seq in 1:length(fasta_results)) {
    for (specie in species_names) {
      current_seq = fasta_results[fasta_seq]
      #print(specie)
      #print(current_seq@ranges@NAMES)
      if (grepl(pattern = specie, x = current_seq@ranges@NAMES, fixed = TRUE)) {
        if(!specie %in% found_species_fasta){
          found_species_fasta <- c(found_species_fasta, specie)
          print(paste("Missing species: ", specie))
        }
      }
    }
  }
  print("Found sequences for: ")
  print(found_species_fasta)
  return(fasta_results)
}

extract_taxonomy_from_taxtable <- function(tax_table, species_names){
  #empty df
  extracted_tax_table <- data.frame(row.names = colnames(colnames(tax_table)))
  
  for (row in 1:nrow(tax_table)) {
    for(specie in species){
      if(tax_table[row,]$species == specie)
        #if(grepl(pattern = specie, x = taxonomy[row,]$species))  
        extracted_tax_table <- rbind(extracted_tax_table, tax_table[row,])
    }
  }
  
  extracted_tax_table <- extracted_tax_table[order(extracted_tax_table$species, decreasing = FALSE),]
  
  print("Sequences extraction complete")
  # Check if all species were found in taxonomy
  found_species_tax <- c()
  
  for (tax_spp in 1:nrow(extracted_tax_table)) {
    for (specie in species) {
      if (specie == extracted_tax_table[tax_spp,]$species) {
        if(!specie %in% found_species_tax){
          found_species_tax <- c(found_species_tax, specie)
          print(paste("Missing species: ", specie))
        }
      }
    }
  }
  print("Found sequences for: ")
  print(found_species_fasta)
}

# Species names for searching in soil microbiome species in fasta file DB.
# Bacillus megaterium is now Priestia megaterium, in taxonomy table name is correct.
species = c("Acidovorax_delafieldii", "Arthrobacter_humicola", "Bacillus_altitudinis", "Bacillus_subtilis",
            "Flavobacterium_pectinovorum", "Bacillus_megaterium", "Pseudomonas_koreensis", "Rhodopseudomonas_palustris")

soil_extract_fasta <- extract_seqs_from_fasta_name(fasta_file = data_fasta, species_names = species)

# Save fasta file
writeXStringSet(x = fasta_results, filepath = "C:/Users/Desktop/species_taxid.fasta")

### Extracting from taxonomy table
# Species name for extracting from taxonomy table
species = c("Acidovorax delafieldii", "Arthrobacter humicola", "Bacillus altitudinis", "Bacillus subtilis",
            "Flavobacterium pectinovorum", "Priestia megaterium", "Pseudomonas koreensis", "Rhodopseudomonas palustris")

soil_extract_taxtable <- extract_taxonomy_from_taxtable(tax_table = data_fasta, species_names = species)

write.table(soil_extract_taxtable, file = "C:/Users/Desktop/taxonomy.tsv")
