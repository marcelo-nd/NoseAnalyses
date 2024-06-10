# Upzip Sequences

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

### rRNA

# Unzip sequences
export SC100_Scree=/mnt/f/SequencingData/SynCom100

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC100_Scree/Screening -o $SC100_Scree

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC100_Scree/fastq -o $SC100_Scree

# emu

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/RRN_db

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC100_Scree/fastq_qc -o $SC100_Scree -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa16copies.csv

Rscript $EMUWRAPPER_LOC/tables_merging.R $Srrna/emu_results "TRUE" /mnt/f/16Sdatabases/LaCa16copies.csv
