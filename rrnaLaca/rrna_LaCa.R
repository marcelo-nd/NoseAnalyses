#install.packages("renv")
#library("renv")

#renv::init("C:/Users/marce/Desktop/rrnaLaca")

renv::install("metacoder")
renv::install("stringr")

renv::snapshot()

renv::restore()

tst_fasta <- metacoder::read_fasta("C:/Users/marce/Desktop/C_pro.fasta")


pcr_result <- metacoder::primersearch_raw(input = tst_fasta,
                                          forward = c("rrnaF" = "AGRGTTYGATYHTGGCTCAG"),
                                          reverse = c("rrnaR" = "ACCRCCCCAGTHRAACT"),
                                          mismatch = 10)

nchar(pcr_result$amplicon[4])

# Get all genomes
# Define the directory path
directory_path <- "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Nasal Genomes/LaCaRenamed"

# Get the list of files with .fasta extension
fasta_files <- list.files(path = directory_path, pattern = "\\.fasta$", full.names = TRUE)

rrnaF = "AGRGTTYGATYHTGGCTCAG"
rrnaR = "ACCRCCCCAGTHRAACT"

for (fasta_file in fasta_files[1:10]) {
  # Print the list of .fasta files
  print(fasta_file)
  # Get pcr result
  # Get amplicons/filter by size
  # add to sequences fasta file
}

pcr_result2 <- metacoder::primersearch_raw(input = metacoder::read_fasta(fasta_files[10]),
                                          forward = rrnaF,
                                          reverse = rrnaR,
                                          mismatch = 10)
nchar(pcr_result2$amplicon[13])

pcr_result2$amplicon[1]

