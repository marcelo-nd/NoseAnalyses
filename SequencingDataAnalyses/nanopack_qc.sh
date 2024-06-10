conda create --name=nanopack python

conda activate nanopack

pip install NanoPlot

# or pip install NanoPlot --upgrade

pip install NanoComp

pip install nanoQC

mkdir chopper0.8

cd chopper0.8

wget https://github.com/wdecoster/chopper/releases/download/v0.8.0/chopper-linux.zip

unzip chopper-linux.zip

rm chopper-linux.zip

chmod +x chopper

# add chopper to path
if [ -d "$HOME/chopper0.8" ] ; then
  PATH="$PATH:$HOME/chopper0.8"
fi

NanoPlot --fastq /mnt/f/SequencingData/Tests/1_16S_First_Test_170124/no_sample/20240117_1510_MN45148_FAY00412_e8d22899/fastq_pass/barcode01/fastq/barcode01_concat.fastq -o /mnt/f/SequencingData/Tests/1_16S_First_Test_170124/ 

/mnt/f/SequencingData/ReproData/1_16S_First_Test_170124_repro/basecalling/pass/barcode01

# Unzip sequences

sed -i 's/\r//' /mnt/c/Users/marce/Documents/GitHub/EmuWrapper/emu_wrapper_unzipper.sh

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export SC1=/mnt/f/SequencingData/ReproData/1_16S_First_Test_170124_repro/basecalling/pass

. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $SC1 -o $SC1/fastaq

####
NanoPlot --fastq /mnt/f/SequencingData/ReproData/1_16S_First_Test_170124_repro/basecalling/pass/barcode01/fastq/barcode01_concat.fastq -o /mnt/f/SequencingData/ReproData/1_16S_First_Test_170124_repro/

############
# QC
sed -i 's/\r//' /mnt/c/Users/marce/Documents/GitHub/EmuWrapper/emu_wrapper_qc.sh

EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper"

export SC1=/mnt/f/SequencingData/ReproData/8_confirmation_rrna_repro_tests/basecalling/pass

. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $SC1 -o $SC1/fastq_qc

####
# Barcode 8
NanoPlot --fastq /mnt/f/SequencingData/Tests/8_confirmation_rrna_230324/no_sample/20240323_1548_MN45148_FAY35039_f0c3d280/fastq_pass/barcode01/fastq/barcode01_concat.fastq -o /mnt/f/SequencingData/Tests/8_confirmation_rrna_230324

NanoPlot --fastq /mnt/f/SequencingData/ReproData/8_confirmation_rrna_repro/basecalling/pass/barcode01/fastq/barcode01_concat.fastq -o /mnt/f/SequencingData/ReproData/8_confirmation_rrna_repro/

NanoPlot --fastq /mnt/f/SequencingData/ReproData/8_confirmation_rrna_repro_tests/basecalling/pass/barcode01/fastq_qc/barcode01_qc.fastq -o /mnt/f/SequencingData/ReproData/8_confirmation_rrna_repro_tests/qc_nanoplot

# rrna tests

NanoPlot --fastq /mnt/f/SequencingData/NanoporeTech/rrnaTests/fastq/barcode01/barcode01_concat.fastq -o /mnt/f/SequencingData/NanoporeTech/rrnaTests/

NanoPlot --fastq /mnt/f/SequencingData/NanoporeTech/rrnaTests/fastq_qc/barcode01/barcode01_qc.fastq -o /mnt/f/SequencingData/NanoporeTech/rrnaTests/

NanoPlot --fastq /mnt/f/SequencingData/NanoporeTech/rrnaTests/fastq_qc2/barcode01/barcode01_qc.fastq -o /mnt/f/SequencingData/NanoporeTech/rrnaTests/

NanoPlot --fastq /mnt/f/SequencingData/NanoporeTech/rrnaTests/fastq_qc_only_short/barcode01/barcode01_qc.fastq -o /mnt/f/SequencingData/NanoporeTech/rrnaTests/

