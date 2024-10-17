# Upzip Sequences

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

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
export S27FII=/mnt/f/SequencingData/NanoporeTech/27FII

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $S27FII -o $S27FII

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $S27FII/fastq -o $S27FII -q 10 -l 800 -h 2500

# emu
# ZCS and Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq_qc -o $S27FII -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03,barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv

# Soil SynCom with Karo DB
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/Karo_plus_050924

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq_qc -o $S27FII -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "FALSE" -p /mnt/f/16Sdatabases/SoilSC_copies.csv

### V34
# Unzip sequences
export V34=/mnt/f/SequencingData/NanoporeTech/V34

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $V34 -o $V34

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $V34/fastq -o $V34 -q 10 -l 100 -h 400

# emu
# ZCS and Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq_qc -o $V34 -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03,barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv

# Soil SynCom with Karo DB
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/Karo_plus_050924

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq_qc -o $V34 -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "FALSE" -p /mnt/f/16Sdatabases/SoilSC_copies.csv


### rRNA
# Unzip sequences
export Srrna=/mnt/f/SequencingData/NanoporeTech/rrna

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $Srrna -o $Srrna

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $Srrna/fastq -o $Srrna -q 10 -l 1000 -h 5000

# emu
# ZCS and Nasal SynCom with LaCa DB
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -b "barcode01,barcode02,barcode03,barcode04,barcode05,barcode06" -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv

# Soil SynCom with Karo DB
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/Karo_plus_050924

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -b "barcode07,barcode08,barcode09" -c "FALSE" -p /mnt/f/16Sdatabases/SoilSC_copies.csv
