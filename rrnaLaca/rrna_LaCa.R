#install.packages("renv")
#library("renv")

#renv::init("C:/Users/marce/Desktop/rrnaLaca")

renv::install("metacoder")
renv::install("stringr")

renv::snapshot()

renv::restore()

# Get all genomes
# Define the directory path
# LaCa Genomes
directory_path <- "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Nasal Genomes/LaCaRenamed"

# Zymo Community Standard genomes
directory_path <- "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Nasal Genomes/ZymoBIOMICS.STD.refseq.v2/Genomes"

# Remaining Genomes
directory_path <- "C:/Users/marce/Desktop/remGenomes"

# Get the list of files with .fasta extension
fasta_files <- list.files(path = directory_path, pattern = "\\.fasta$", full.names = TRUE)

rrnaF = "AGRGTTYGATYHTGGCTCAG"
rrnaR = "ACCRCCCCAGTHRAACT"

amplicon_fasta_list <- c()

for (fasta_file in fasta_files) {
  # Print the current fasta file
  #print(fasta_file)
  # Split the string by "/"
  parts <- strsplit(fasta_file, "/")[[1]]
  # Take the last part of the split string
  last_part <- tail(parts, n = 1)
  # Remove ".fasta" from the last part
  bacteria <- sub("\\.fasta$", "", last_part)
  
  # Print the result
  print(bacteria)
  
  # Get pcr result
  pcr_result <- metacoder::primersearch_raw(input = metacoder::read_fasta(fasta_file),
                                            forward = rrnaF,
                                            reverse = rrnaR,
                                            mismatch = 10)
  # print number of amplicons
  #print(nrow(pcr_result))
  # Get amplicons/filter by size. FIltering according to MIrROR paper, 3500 to 7000 bp.
  amplicons <- pcr_result[nchar(pcr_result[["amplicon"]]) >= 3500 & nchar(pcr_result[["amplicon"]]) <= 7000, ]$amplicon
  #print(length(amplicons))
  # add to sequences fasta file
  if (length(amplicons) > 0) {
    amplicon_counter <- 1
    for (amplicon in amplicons) {
      
      # check if amplicon contains unallowed characters
      allowed_chars <- "ACGT"
      
      # Create a regular expression pattern for allowed characters
      pattern <- paste0("[^", allowed_chars, "]")
      
      # Check if the string contains any character not in the allowed subset
      contains_disallowed <- grepl(pattern, amplicon)
      
      # Print the result
      if (contains_disallowed) {
        cat("The string contains disallowed characters!!!!!!!!!!!!!!!\n")
      } else {
        cat("The string only contains allowed characters.\n")
      }
      #### end of checking
      
      amplicon_fasta_list <- c(amplicon_fasta_list, paste0(">", paste(bacteria, amplicon_counter, sep = "_")))
      amplicon_fasta_list <- c(amplicon_fasta_list,amplicon)
      amplicon_counter <- amplicon_counter + 1
    }
  }
}

print(amplicon_fasta_list)

writeLines(amplicon_fasta_list, con = "C:/Users/marce/Desktop/Zymo_rRNA.fasta")

