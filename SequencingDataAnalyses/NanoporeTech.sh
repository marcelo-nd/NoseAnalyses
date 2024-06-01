# Upzip Sequences

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

### 27F

export S27F=/mnt/f/SequencingData/NanoporeTech/27FTests

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $S27F -o $S27F

### 27F-II

### V34

### rRNA

# Unzip sequences
export Srrna=/mnt/f/SequencingData/NanoporeTech/rrnaTests

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $Srrna -o $Srrna

# QC

export Srrna=/mnt/f/SequencingData/NanoporeTech/rrnaTests

. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $Srrna/fastq -o $Srrna

# emu

export Srrna=/mnt/f/SequencingData/NanoporeTech/rrnaTests

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $Srrna/fastq_qc2 -o $Srrna -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv
