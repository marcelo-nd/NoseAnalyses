# create conda environment for emu
conda create --name emu_py37 python=3.7

# install emu
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install emu

# R steps
# Install R
sudo apt install r-base-core

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

# SynCom Test 5 (05.03.24)
conda activate emu_py37

export SC5=/mnt/e/1_NoseSynComProject/SequencingData/5_240305_Primer_27F-Test

. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "FALSE" -s $SC5/no_sample/20240305_1407_MN45148_FAY00412_2df0c6ae/fastq_pass -o $SC5/results -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv

# SynCom Test 6
conda activate emu_py37

export SC6=/mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna

. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "TRUE" -s $SC6/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass -o $SC5/results -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv


./MIrROR.py -V -d DBDIR -o ./result /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode06/fastq/barcode06_concat.fastq

# SynCom Test 7 (11.03.24)
conda activate emu_py37

export SC7=/mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks
. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "TRUE" -s $SC7/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/fastq_pass -o $SC7/results  -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv


# SynCom Test 7 (11.03.24) (part 2, long reads only decompressing)
conda activate emu_py37

export SC7=/mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks
. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "TRUE" -s $SC7/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2 -o $SC7/results_long  -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv


# SynCom Test 8 (03.04.24) (16s kit test of Mock communities)
conda activate emu_py37

export SC9=/mnt/e/1_NoseSynComProject/SequencingData/9_16STest_MockDNA_030424
. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "FALSE" -s $SC9/no_sample/20240403_1526_MN45148_FAY00412_a67f721c/fastq_pass -o $SC9/results_long  -c "TRUE" -p /mnt/c/LaCa16copies.csv
