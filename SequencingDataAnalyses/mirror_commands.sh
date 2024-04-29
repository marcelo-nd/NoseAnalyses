EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

MIRRORWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/MIrROR_Wrapper"

##### First experiment

~/MIrROR/MIrROR.py -V -d DBDIR -o ./result8 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode06/fastq/barcode06_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result9 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode07/fastq/barcode07_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result10 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode08/fastq/barcode08_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result11 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode09/fastq/barcode09_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result12 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode10/fastq/barcode10_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result13 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode11/fastq/barcode11_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result7 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode12/fastq/barcode12_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result14 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/6_070324_Creat_rrna/no_sample/20240307_1134_MN45148_FAY35039_5bd01a44/fastq_pass/barcode13/fastq/barcode13_concat.fastq


export SC6=/mnt/e/1_NoseSynComProject/SequencingData/6_Creat_rrna_070324/results_mirror

Rscript $MIRRORWRAPPER_LOC/tables_merging.R $SC6 "TRUE" /mnt/c/LaCa16copies.csv

##### Second experiment
export SC7=/mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/results_mirror

./MIrROR.py -V -d DBDIR -o ./result15 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode11/fastq/barcode11_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result16 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode12/fastq/barcode12_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result17 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode13/fastq/barcode13_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result18 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode14/fastq/barcode14_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result19 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode15/fastq/barcode15_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result20 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode16/fastq/barcode16_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result21 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode17/fastq/barcode17_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result22 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode18/fastq/barcode18_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result23 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode19/fastq/barcode19_concat.fastq

./MIrROR.py -V -d DBDIR -o ./result24 -m 1000 -b 2000 /mnt/e/1_NoseSynComProject/SequencingData/7_PrimerConfTest_Mocks/no_sample/20240311_1213_MN45148_FAY35039_b3637fda/2/barcode20/fastq/barcode20_concat.fastq

# Merge Tables
Rscript $MIRRORWRAPPER_LOC/tables_merging.R $SC7 "TRUE" /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv


##### Third experiment

export SC8=/mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna

export SC8=/mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/results_mirror

# Unzipping, in the future Mirror_wrapper will do this
. $EMUWRAPPER_LOC/emu_wrapper.sh -d $EMU_DATABASE_DIR -z "TRUE" -s $SC8/no_sample/20240323_1548_MN45148_FAY35039_f0c3d280/fastq_pass -o $SC8/results  -c "TRUE" -p /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv

cd /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/no_sample/20240323_1548_MN45148_FAY35039_f0c3d280/fastq_pass

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/result01 -m 1000 -b 2000 ./barcode01/fastq/barcode01_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/result02 -m 1000 -b 2000 ./barcode02/fastq/barcode02_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/result03 -m 1000 -b 2000 ./barcode03/fastq/barcode03_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/result04 -m 1000 -b 2000 ./barcode04/fastq/barcode04_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/result05 -m 1000 -b 2000 ./barcode05/fastq/barcode05_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/results_mirror/result06 -m 1000 -b 2000 ./barcode06/fastq/barcode06_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/results_mirror/result07 -m 1000 -b 2000 ./barcode07/fastq/barcode07_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/results_mirror/result08 -m 1000 -b 2000 ./barcode08/fastq/barcode08_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/results_mirror/result09 -m 1000 -b 2000 ./barcode09/fastq/barcode09_concat.fastq

~/MIrROR/MIrROR.py -V -d ~/MIrROR/DBDIR/ -o /mnt/e/1_NoseSynComProject/SequencingData/8_confirmation_rrna/results_mirror/result10 -m 1000 -b 2000 ./barcode10/fastq/barcode10_concat.fastq

# Merging tables
Rscript $MIRRORWRAPPER_LOC/tables_merging.R $SC8 "TRUE" /mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/SequencingDataAnalyses/LaCa16copies.csv
