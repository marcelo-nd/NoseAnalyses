# create conda environment for emu (chose a new name!)
conda create --name emu python=3.7

conda activate emu

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

# Download Database (instructions from emu GitHub repo)
pip install osfclient
mkdir /home/marcelo/emuDB
export EMU_DATABASE_DIR=/home/marcelo/emuDB #change for your own path
cd ${EMU_DATABASE_DIR}
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
tar -xvf emu.tar

# Set EmuWrapper PATH
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper" # change for your path

#fix line endings of EmuWrapper.sh (this prevents some posterior errors)
sed -i 's/\r//' $EMUWRAPPER_LOC/emu_wrapper.sh

# Now lets install Nanopack chopper package to do the QC
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

### Now we can analyse sequences!

########## Sequencing Run from Karo

# Activate emu if it is not activated
conda activate emu

# Set EmuWrapper PATH
EMUWRAPPER_LOC="/mnt/c/Users/marce/Documents/GitHub/EmuWrapper" # change for your path

# Store sequences path, change for your own PATH
export seqs_paths=/mnt/c/Users/marce/Desktop/pass

# Unzip sequences
. $EMUWRAPPER_LOC/emu_wrapper_unzipper.sh -s $seqs_paths -o $seqs_paths

# Do the QC
. $EMUWRAPPER_LOC/emu_wrapper_qc.sh -s $seqs_paths/fastq -o $seqs_paths -q 10 -l 1000 -h 5000

# Run emu!
export EMU_DATABASE_DIR=/home/marcelo/emuDB # This is the emu Database, change for your own path
# or use Database from Karo's genomes:
export EMU_DATABASE_DIR=/mnt/e/16Sdatabases/Karo_plus_050924 # change for your own path

# Copy number database
export COPY_No_DB=/mnt/e/16Sdatabases/SoilSC_copies.csv # This file contains the data about number of 16S copies in Karo's genomes, change for your own path

. $EMUWRAPPER_LOC/emu_wrapper_run_emu.sh -s $seqs_paths/fastq_qc -o $seqs_paths -d $EMU_DATABASE_DIR -c "TRUE" -p $COPY_No_DB
