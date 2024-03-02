# create conda environment for emu
conda create --name emu_py37 python=3.7

# install emu
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install emu

# R steps
# Install R

# Activate R
R

# Install required packages
install.packages("readr", lib = "/usr/lib/R/library")
install.packages("dplyr", lib = "/usr/lib/R/library")
install.packages("plyr", lib = "/usr/lib/R/library")

# quit R
# quit()


# Create custom database
cd /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database

conda activate emu_py37

emu build-database LaCaEmu --sequences /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/DB_files/LaCa_sequences.fasta --seq2tax /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/DB_files/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/DB_files/LaCa_taxonomy.tsv

sudo apt install r-base-core

# Set EmuWrapper
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export EMU_DATABASE_DIR=./emuDB

export EMU_DATABASE_DIR=/mnt/e/1_NoseSynComProject/SequencingData/LaCaEmu

#fix line endings of EmuWrapper.sh
sed -i 's/\r//' /mnt/c/Users/marce/Documents/GitHub/EmuWrapper/emu_wrapper.sh

# SynCom Test 1 (17.01.24)
conda activate emu_py37

export SC1=/mnt/e/1_NoseSynComProject/SequencingData/1_16S_First_Test_17_01_24

. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "FALSE" -s $SC1/no_sample/20240117_1510_MN45148_FAY00412_e8d22899/fastq_pass -o $SC1/results  -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv

# SynCom Test 2 (24.01.24)
conda activate emu_py37

. emu_wrapper.sh -d /home/marcelo/LaCa -s /mnt/e/NoseSynComProject/SequencingData/2_Test2_240124/no_sample/20240124_1634_MN45148_FAY00412_07f41999/fastq_pass -o /mnt/e/NoseSynComProject/SequencingData/2_Test2_240124/results

# SynCom Test 3 (12.02.24)
conda activate emu_py37

export SC3=/mnt/e/NoseSynComProject/SequencingData/3_SynComTest3

. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "TRUE" -s $SC3/no_sample/20240212_1656_MN45148_FAY00412_e1705dac/fastq_pass -o $SC3/results -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv

# SynCom Test 4 (28.02.24)
conda activate emu_py37

export SC4=/mnt/e/1_NoseSynComProject/SequencingData/Test4-SynComNewPrimers

. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "TRUE" -s $SC4/no_sample/20240228_1510_MN45148_FAY35039_711f6cab/fastq_pass -o $SC4/results -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv
