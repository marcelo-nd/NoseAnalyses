conda activate emu

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

### Screening

# Unzip sequences
export SC100_Scree=/mnt/f/SequencingData/SynCom100/Screening

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC100_Scree/Sequences -o $SC100_Scree

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC100_Scree/fastq -o $SC100_Scree -q 10 -l 1000 -h 5000

# emu

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC100_Scree/fastq_qc -o $SC100_Scree -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv



### Timepoints

# Unzip sequences
export SC100_TP=/mnt/f/SequencingData/SynCom100/TheChampions

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC100_TP -o $SC100_TP

# QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC100_TP/fastq -o $SC100_TP -q 10 -l 1000 -h 5000

# emu

export EMU_DATABASE_DIR=/mnt/f/16Sdatabases/LaCa_rRNA_Emu_210824

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $SC100_TP/fastq_qc -o $SC100_TP -d $EMU_DATABASE_DIR -c "TRUE" -p /mnt/f/16Sdatabases/LaCa_copies.csv
