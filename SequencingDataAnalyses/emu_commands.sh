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


















# Set EmuWrapper
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"



#### SequencingRunBatch5 25.10.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SequencingRunBatch5/no_sample_id/20241025_1411_MN45148_FAZ32262_7e5d6e09/fastq_pass

#### SequencingRunBatch5 28.10.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_SR2/no_sample_id/20241028_1159_MN45148_AVN094_845724a9/fastq_pass

#### SequencingRunBatch5 29.10.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SeqRun3_291024/no_sample_id/20241029_1750_MN45148_AUN642_9d368fc1/fastq_pass

#### SequencingRunBatch5 12.11.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SeqRun_SC100_121124/no_sample_id/20241112_1515_MN45148_auh943_d79d3d47/fastq_pass

#### SequencingRunBatch5 16.11.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_SeqRun_151124/no_sample_id/20241115_1500_MN45148_FBA32257_7de35303/fastq_pass

#### SequencingRunBatch5 28.11.24

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_281124/no_sample_id/20241128_1507_MN45148_FAY35296_ff71f5b3/fastq_pass

#### SequencingRunBatch5 20.01.25

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_200125/no_sample_id/20250120_1209_MN45148_aui012_771ff743/fastq_pass

#### SequencingRunBatch5 24.01.25

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_230125/no_sample_id/20250123_1637_MN45148_AVM817_037deab6/fastq_pass

#### SequencingRunBatch5 30.01.25

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_290125/no_sample_id/20250129_1629_MN45148_AUI367_bf1e5858/fastq_pass

### Duda SC12

export seqs_paths=/mnt/e/SequencingData/SynCom100/duda/SC12

### Duda SC13
export seqs_paths=/mnt/e/SequencingData/SynCom100/duda/SC13

#### SequencingRunBatch5 05.02.25
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_050224/no_sample_id/20250205_1052_MN45148_FBA33889_f4d732c7/fastq_pass

#### SequencingRun 13.02.25

export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_130225/no_sample_id/20250213_1012_MN45148_AYU138_31d9206a/fastq_pass

#### Check 
export seqs_paths=/mnt/e/SequencingData/SynCom100/duda

#### SequencingRun 18.02.25 (nanopore computer)
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_180225/no_sample_id/20250218_1259_MN45148_AWS287_8c3aa7bd/fastq_pass


#### SequencingRun 31.01.25 (nanopore computer)
EMUWRAPPER_LOC="/mnt/c/Users/AGLinkNanopore/Documents/Marcelo/EmuWrapper"
export seqs_paths=/mnt/c/data/Marcelo/SC100_300125/no_sample_id/20250131_1506_MN45148_awo447_fd9b78e9/fastq_pass

#### SequencingRun 07.02.25 (nanopore computer)
EMUWRAPPER_LOC="/mnt/c/Users/AGLinkNanopore/Documents/Marcelo/EmuWrapper"
export seqs_paths=/mnt/c/data/Marcelo/SC100_070225/no_sample_id/20250207_1300_MN45148_FAY35296_64b2d4a6/fastq_pass

#### SequencingRun 14.02.25 (nanopore computer)
EMUWRAPPER_LOC="/mnt/c/Users/AGLinkNanopore/Documents/Marcelo/EmuWrapper"
export seqs_paths=/mnt/c/data/Marcelo/SC100_140225/no_sample_id/20250214_1618_MN45148_AYY707_b38a097e/fastq_pass

#### SequencingRun 21.02.25 (nanopore computer)
EMUWRAPPER_LOC="/mnt/c/Users/AGLinkNanopore/Documents/Marcelo/EmuWrapper"
export seqs_paths=/mnt/c/data/Marcelo/SC100_200225/no_sample_id/20250220_1746_MN45148_FAZ32262_c1bd7e7d/fastq_pass

#### SequencingRun 21.02.25 (nanopore computer)
EMUWRAPPER_LOC="/mnt/c/Users/AGLinkNanopore/Documents/Marcelo/EmuWrapper"
export seqs_paths=/mnt/c/data/Marcelo/SC100_200225/no_sample_id/20250220_1746_MN45148_FAZ32262_c1bd7e7d/fastq_pass

#### SequencingRun 07.03.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/NS_SC100_070325/no_sample_id/20250307_1222_MN45148_FBA33889_f5fe9dab/fastq_pass

#### SequencingRun 08.04.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SCPlus_070425/no_sample_id/20250407_1236_MN45148_AYT875_52253316/fastq_pass

#### SequencingRun 17.04.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100Plus_try2_160425/no_sample_id/20250416_1346_MN45148_ayz085_541b5a66/fastq_pass

#### SequencingRun 25.04.25
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_PaoloSC_250425/no_sample_id/20250425_1330_MN45148_AWM177_48528f9c/fastq_pass

# Emu for run 25.04.25
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 300 -h 5000

export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_rRNA_Emu_181124
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03,barcode04,barcode05,barcode06,barcode07,barcode08,barcode09,barcode10,barcode23" -c "TRUE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_copies.csv

export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC_mod
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -b "barcode11,barcode12,barcode13,barcode14,barcode15,barcode16,barcode17,barcode18,barcode19,barcode20,barcode21,barcode22" -c "TRUE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC/Ara_Library_copies.csv

#### Combined run 17.04.25 and 25.04.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/Comb_SC100_PaoloSC_250425_KemenSC_260425
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 300 -h 5000
export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC_mod
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -b "barcode11,barcode12,barcode13,barcode14,barcode15,barcode16,barcode17,barcode18,barcode19,barcode20,barcode21,barcode22" -c "TRUE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC/Ara_Library_copies.csv


#### Run 15.05.25 (SEVERAL PRIMER FROM BARCODES 1-9)
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/Paolo_SC100_140525/no_sample_id/20250514_1614_MN45148_AWR003_8e853f3f/fastq_pass
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 150 -h 1100
export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC_mod
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03,barcode04,barcode05,barcode06,barcode07,barcode08,barcode09" -c "TRUE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC/Ara_Library_copies.csv


#### Run Comb_paolo 15.05.25 (SEVERAL PRIMER FROM BARCODES 1-9)
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/Comb_Paolo_140525_150525
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 150 -h 1100
export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC_mod
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03,barcode04,barcode05,barcode06,barcode07,barcode08,barcode09" -c "FALSE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/KemenSC/Ara_Library_copies.csv


#### Run SC100 15.05.25 (SEVERAL PRIMER FROM BARCODES 1-9)
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/Comb_SC100_140525_150525
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 1000 -h 5000
export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_rRNA_Emu_181124
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR  -c "TRUE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_copies.csv

"barcode12,barcode13,barcode14,barcode15,barcode16,barcode17,barcode18,barcode19,barcode20,barcode21,barcode22" 

#### Run NPT 23.05.25
EMUWRAPPER_LOC="/mnt/c/Users/AGLinkNanopore/Documents/Marcelo/EmuWrapper"
export seqs_paths=/mnt/c/data/Marcelo/NanoporeTech_230525/no_sample_id/20250523_1713_MN45148_FBB06794_1bd5bf2c/fastq_pass





EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 1000 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_rRNA_Emu_181124
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "FALSE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_copies.csv

# emu (nanopore computer)
export EMU_DATABASE_DIR=/mnt/c/Users/AGLinkNanopore/Documents/Marcelo/16Sdatabases/LaCa_rRNA_Emu_181124
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "FALSE" -p /mnt/c/Users/AGLinkNanopore/Documents/Marcelo/16Sdatabases/LaCa_copies.csv


# plant database
emu build-database KemenSC_DB --sequences /mnt/d/2_OtherProjects/Paolo_Sequencing/DB_files/plant_syncom.fasta --seq2tax /mnt/d/2_OtherProjects/Paolo_Sequencing/DB_files/plant_seq2tax.map --taxonomy-list /mnt/d/2_OtherProjects/Paolo_Sequencing/DB_files/plant_taxonomy.tsv


# Only Laca
emu build-database LaCa_DB --sequences /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/only_LaCa/LaCa_rRNA.fasta --seq2tax /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/only_LaCa/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/only_LaCa/LaCa_taxonomy.tsv

# Only Zymo
emu build-database Zymo_DB --sequences /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/only_zymo/LaCa_rRNA.fasta --seq2tax /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/only_zymo/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/marce/'OneDrive - UT Cloud'/'1_NoseSynCom Project'/'Nasal Genomes'/LaCa_rRNA/only_zymo/LaCa_taxonomy.tsv



#### Run NPT 28.05.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/NT_28052025/no_sample_id/20250528_1532_MN45148_FBA33889_840e0bfa/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 200 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/RRN_db
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "FALSE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_copies.csv



#### Run NPT 28.05.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/d/1_NoseSynComProject/SequencingData/OriginalRuns/SC100-Repetition_Experiment/no_sample_id/20250731_1603_MN45148_FAY35296_2d1229a5/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 1000 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_rRNA_Emu_181124
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_copies.csv


#### Run NPT 04.07.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/Rep_exp_2/no_sample_id/20250804_1145_MN45148_FBA33889_6231beab/basecalling/pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 1000 -h 5000

# emu
export EMU_DATABASE_DIR=/mnt/e/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_rRNA_Emu_181124
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/16Sdatabases/16s_copies.csv


#### Run NPT 20.008.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/Nanoporetech_200825/no_sample_id/20250820_1858_MN45148_FBB06794_bd0ec299/basecalling/pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 500 -h 1100

# emu
export EMU_DATABASE_DIR=/mnt/e/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_rRNA_Emu_181124
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "FALSE" -p /mnt/d/16Sdatabases/16s_copies.csv



#### Run AAs Exp 28.08.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/AAsExp_270825/no_sample_id/20250827_1852_MN45148_FBB06794_a4f73dd9/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 500 -h 5500

# emu
export EMU_DATABASE_DIR=/mnt/e/1_NoseSynComProject/SequencingData/16Sdatabases/LaCa_rRNA_Emu_181124
export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/nose_sc_db_200824
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/16Sdatabases/16s_copies.csv



#### Run Paolo 28.08.25
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/PaoloBatch2_Bact1/no_sample_id/20250912_1552_MN45148_FBB05369_22cef7a6/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 500 -h 1500

# emu
export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/plant_zymo_121625
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/16Sdatabases/16s_copies.csv



#### 
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/PaoloSC2_221025_Nose/no_sample_id/20251022_1402_MN45148_FBA97141_ce19b3fa/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 500 -h 1500

# emu
export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/nose_sc_db_200824
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/16Sdatabases/16s_copies.csv





#### 
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/PaoloSC2_221025_Paolo/no_sample_id/20251022_1402_MN45148_FBA97141_ce19b3fa/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 500 -h 1500

# emu
export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/plant_sc_db_061025

export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/16Sdatabases/16s_copies.csv




#### 
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/PaoloSC_011125/no_sample_id/20251101_1501_MN45148_FBE49262_d86ea796/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 500 -h 1500

# emu
export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/plant_sc_300425

export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/16Sdatabases/16s_copies.csv

# barcode11 barcode17



#### 
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"
export seqs_paths=/mnt/e/1_NoseSynComProject/SequencingData/OriginalRuns/SC100_Cocultures/no_sample_id/20251030_1310_MN45148_FBA97141_354f78d2/fastq_pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 500 -h 1500

# emu
export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/plant_sc_db_061025

export EMU_DATABASE_DIR=/mnt/d/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/d/16Sdatabases/16s_copies.csv