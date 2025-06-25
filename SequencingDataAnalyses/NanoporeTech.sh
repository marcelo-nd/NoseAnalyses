EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

# Genral Database
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/RRN_db

# Upzip Sequences
### 27F
# Unzip sequences
export S27F=/mnt/f/SequencingData/NanoporeTech/27F

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $S27F -o $S27F

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $S27F/fastq -o $S27F

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCaEmu

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27F/fastq_qc -o $S27F -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv



### 27F-II
# Unzip sequences
export S27FII=/mnt/e/SequencingData/NanoporeTech/1_Sequences/27FII

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $S27FII -o $S27FII

# QC (~1492)
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $S27FII/fastq -o $S27FII -q 10 -l 800 -h 1600 2>&1 | tee $S27FII/fastq_qc/qc_output.txt

# emu
# General database
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq_qc -o $S27FII -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq_qc -o $S27FII -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq_qc -o $S27FII -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Gut SynCom20 with Gut SC20 DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/gut_sc_db_260225

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq_qc -o $S27FII -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq_qc -o $S27FII -d $EMU_DATABASE_DIR -b "barcode10,barcode11,barcode12" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv


### V34 (~444)
# Unzip sequences
export V34=/mnt/e/SequencingData/NanoporeTech/1_Sequences/V34

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $V34 -o $V34

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $V34/fastq -o $V34 -q 10 -l 250 -h 500 2>&1 | tee $V34/fastq_qc/qc_output.txt

# emu
# General database
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq_qc -o $V34 -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq_qc -o $V34 -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq_qc -o $V34 -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Gut SynCom20 with GutSC20 DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/gut_sc_db_260225

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq_qc -o $V34 -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq_qc -o $V34 -d $EMU_DATABASE_DIR -b "barcode10,barcode11,barcode12" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv


### rRNA (~7000)
# Unzip sequences
export Srrna=/mnt/e/SequencingData/NanoporeTech/1_Sequences/rrna

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $Srrna -o $Srrna

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $Srrna/fastq -o $Srrna -q 10 -l 2000 -h 7000 2>&1 | tee $Srrna/fastq_qc/qc_output.txt

# emu
# General database
. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Gut SynCom20 with rRNA Mirror DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/gut_sc_db_260225

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -b "barcode10,barcode11,barcode12" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

### v13 (507bp)
# Unzip sequences
export V13=/mnt/e/SequencingData/NanoporeTech/1_Sequences/v13

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $V13 -o $V13

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $V13/fastq -o $V13 -q 10 -l 300 -h 600 2>&1 | tee $V13/fastq_qc/qc_output.txt

# emu
# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V13/fastq_qc -o $V13 -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V13/fastq_qc -o $V13 -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Gut SynCom20 with rRNA Mirror DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/gut_sc_db_260225

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V13/fastq_qc -o $V13 -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V13/fastq_qc -o $V13 -d $EMU_DATABASE_DIR -b "barcode10,barcode11,barcode12" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv


### v38 (~1037 bp)
# Unzip sequences
export V38=/mnt/e/SequencingData/NanoporeTech/1_Sequences/v38

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $V38 -o $V38

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $V38/fastq -o $V38 -q 10 -l 500 -h 1100 2>&1 | tee $V38/fastq_qc/qc_output.txt

# emu
# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V38/fastq_qc -o $V38 -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V38/fastq_qc -o $V38 -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Gut SynCom20 with rRNA Mirror DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/gut_sc_db_260225

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V38/fastq_qc -o $V38 -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V38/fastq_qc -o $V38 -d $EMU_DATABASE_DIR -b "barcode10,barcode11,barcode12" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

### v79 (~377 bp)
# Unzip sequences
export V79=/mnt/e/SequencingData/NanoporeTech/1_Sequences/v79

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $V79 -o $V79

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $V79/fastq -o $V79 -q 10 -l 200 -h 400 2>&1 | tee $V79/fastq_qc/qc_output.txt

# emu
# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V79/fastq_qc -o $V79 -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/nose_sc_db_200824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V79/fastq_qc -o $V79 -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Gut SynCom20 with rRNA Mirror DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/gut_sc_db_260225

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V79/fastq_qc -o $V79 -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V79/fastq_qc -o $V79 -d $EMU_DATABASE_DIR -b "barcode10,barcode11,barcode12" -c "TRUE" -p /mnt/e/16Sdatabases/16s_copies.csv


### ITS12 (~ 320 bp)
# Unzip sequences
export ITS12=/mnt/e/SequencingData/NanoporeTech/1_Sequences/ITS12

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $ITS12 -o $ITS12

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $ITS12/fastq -o $ITS12 -q 10 -l 100 -h 400 2>&1 | tee $ITS12/qc_output.txt

# emu
# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $ITS12/fastq_qc -o $ITS12 -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "FALSE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $ITS12/fastq_qc -o $ITS12 -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "FALSE" -p /mnt/e/16Sdatabases/16s_copies.csv



### ITS14 (~550 bp)
# Unzip sequences
export ITS14=/mnt/e/SequencingData/NanoporeTech/1_Sequences/ITS14

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $ITS14 -o $ITS14

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $ITS14/fastq -o $ITS14 -q 10 -l 300 -h 600 2>&1 | tee $ITS14/qc_output.txt

# emu
# ZCS with ZCS DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/zymo_db_220525

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $ITS14/fastq_qc -o $ITS14 -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03" -c "FALSE" -p /mnt/e/16Sdatabases/16s_copies.csv

# Plant SynCom with KemenSC DB
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/plant_sc_300425

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $ITS14/fastq_qc -o $ITS14 -d $EMU_DATABASE_DIR -b "barcode04,barcode05,barcode06" -c "FALSE" -p /mnt/e/16Sdatabases/16s_copies.csv
