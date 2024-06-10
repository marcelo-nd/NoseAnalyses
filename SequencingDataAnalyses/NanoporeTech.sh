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

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27F/fastq -o $S27F -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

# Table Merging
Rscript $EMUWRAPPER_LOC/tables_merging.R $S27F/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv


### 27F-II
# Unzip sequences
export S27FII=/mnt/f/SequencingData/NanoporeTech/27FII

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $S27FII -o $S27FII

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $S27FII/fastq -o $S27FII

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCaEmu

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $S27FII/fastq -o $S27FII -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

# Table Merging
Rscript $EMUWRAPPER_LOC/tables_merging.R $S27FII/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv


### V34
# Unzip sequences
export V34=/mnt/f/SequencingData/NanoporeTech/V34

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $V34 -o $V34

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $V34/fastq -o $V34

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCaEmu

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $V34/fastq -o $V34 -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

# Table Merging
Rscript $EMUWRAPPER_LOC/tables_merging.R $V34/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv


### rRNA
# Unzip sequences
export Srrna=/mnt/f/SequencingData/NanoporeTech/rrna

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $Srrna -o $Srrna

# QC
export Srrna=/mnt/f/SequencingData/NanoporeTech/rrnaTests

. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $Srrna/fastq -o $Srrna

# emu
export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

# Table Merging
Rscript $EMUWRAPPER_LOC/tables_merging.R $Srrna/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv












##############################################
### rRNA tests

# Unzip sequences
export Srrna=/mnt/f/SequencingData/NanoporeTech/rrnaTests

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $Srrna -o $Srrna

# QC
export Srrna=/mnt/f/SequencingData/NanoporeTech/rrnaTests

. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $Srrna/fastq -o $Srrna

# emu

export Srrna=/mnt/f/SequencingData/NanoporeTech/rrnaTests

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_medium -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_long -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc_only_short -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

Rscript $EMUWRAPPER_LOC/tables_merging.R $Srrna/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv


########################################## SynCom Batch 1+2

# Unzip sequences
export scb1_2=/mnt/f/SequencingData/SynComExp070624/no_sample/20240607_1513_MN45148_FAY35296_77d43522/fastq_pass

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $scb1_2 -o /mnt/f/SequencingData/SynComExp070624/

# QC
export scb1_2=/mnt/f/SequencingData/SynComExp070624

. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $scb1_2/fastq -o $scb1_2

# Nanoplot

NanoPlot --fastq /mnt/f/SequencingData/SynComExp070624/fastq_qc/barcode01/barcode01_qc.fastq -o /mnt/f/SequencingData/SynComExp070624/

# emu

export scb1_2=/mnt/f/SequencingData/SynComExp070624

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $scb1_2/fastq_qc -o $scb1_2 -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv
