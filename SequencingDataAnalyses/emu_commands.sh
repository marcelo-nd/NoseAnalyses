# create conda environment for emu
conda create --name emu python=3.7

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
install.packages("readr")
install.packages("dplyr")
install.packages("plyr")

# quit R
quit()


# Create custom database
cd /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database

conda activate emu

emu build-database LaCaEmu --sequences /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/DB_files/LaCa_sequences.fasta --seq2tax /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/DB_files/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/DB_files/LaCa_taxonomy.tsv

emu build-database LaCa_rRNA_Emu --sequences /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/DB_files/LaCa_rRNA.fasta --seq2tax /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/DB_files/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/DB_files/LaCa_taxonomy.tsv

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


# SynCom Test 8 (03.04.24) (16s kit test of Mock communities)
conda activate emu_py37

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

export SC10=/mnt/f/SequencingData/OriginalRuns/ReproData/SynComTFBatch1and2_250524

# Unzip sequences

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC10/Basecalling -o $SC10

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC10/fastq -o $SC10

# emu

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC10/fastq_qc -o $SC10 -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv








sed -i 's/\r//' $EMUWRAPPER_LOC/emu_wrapper.sh

. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "TRUE" -s $SC10/no_sample/20240525_1628_MN45148_FAY35039_2aa0c353/fastq_pass -o $SC10/results_long  -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

#################### SC Karo 16S (25.07.24)

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

export SC10=/mnt/f/SequencingData/NanoporeTech/SC_soil

# Unzip sequences

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC10/basecalling/pass -o $SC10

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC10/fastq -o $SC10

# emu

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC10/fastq_qc -o $SC10 -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

Rscript $EMUWRAPPER_LOC/tables_merging.R $SC10/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv

########## Karo multiple barcodes

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export EMU_DATABASE_DIR=/mnt/d/SequencingTemp/16Sdatabases/RRN_db

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/Karo_Emu_200824

export seqs_paths=/mnt/f/SequencingData/Karo/

# Unzip sequences

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC10 -o $SC10

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC10/fastq -o $SC10

# emu

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "FALSE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

Rscript $EMUWRAPPER_LOC/tables_merging.R $SC10/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv


########## Nose tests

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export EMU_DATABASE_DIR=/mnt/d/SequencingTemp/16Sdatabases/LaCa_rRNA

export SC10=/mnt/d/SequencingTemp/KaroSC_NasalSC/no_sample/20240726_1242_MN45148_aun678_1f2572b2/Nose

# Unzip sequences

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC10 -o $SC10

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC10/fastq -o $SC10

# emu

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC10/fastq_qc -o $SC10 -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/SequencingTemp/16Sdatabases/LaCa16copies.csv

Rscript $EMUWRAPPER_LOC/tables_merging.R $SC10/emu_results "TRUE" /mnt/d/SequencingTemp/16Sdatabases/LaCa16copies.csv


#### Karo2

########## Karo multiple barcodes

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export EMU_DATABASE_DIR=/mnt/d/SequencingTemp/16Sdatabases/RRN_db

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/Karo_Emu_200824

export seqs_paths=/mnt/f/SequencingData/Karo2/

# Unzip sequences

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths

# emu

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "FALSE" -p /mnt/f/16Sdatabases/LaCa16copies.csv


# Karos test with DB plus
export SEQS=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/20240829_Karo_PRimer_SoilSynCom2/no_sample/20240829_1253_MN45148_AVM806_08a69708/fastq_pass
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SEQS/fastq_qc -b "barcode07,barcode08,barcode09,barcode10,barcode11,barcode12,barcode13,barcode14,barcode15,barcode16,barcode17,barcode18" -o $SEQS -d $EMU_DATABASE_DIR -c "FALSE"

# Screenind and missing samples

export seqs_paths=/mnt/f/NasalSC100_300924/no_sample_id/20240930_1731_MN45148_FBA32257_3c181886/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 200 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv


#### SequencingRunBatch5 25.10.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SequencingRunBatch5/no_sample_id/20241025_1411_MN45148_FAZ32262_7e5d6e09/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 200 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths/nose -b "barcode01,barcode02,barcode03,barcode04,barcode05,barcode06,barcode07,barcode08,barcode09,barcode10,barcode11,barcode12,barcode13,barcode14,barcode24" -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths/gsc -b "barcode15,barcode16,barcode17,barcode18,barcode19,barcode20,barcode21,barcode22,barcode23" -d $EMU_DATABASE_DIR -c "FALSE" -p /mnt/f/16Sdatabases/LaCa_copies.csv


#### SequencingRunBatch5 28.10.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_SR2/no_sample_id/20241028_1159_MN45148_AVN094_845724a9/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 1000 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv


#### SequencingRunBatch5 29.10.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/Temp/SeqRun3_291024/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 1000 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv
