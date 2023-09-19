#!/bin/bash
echo "Database path: $1"
# ./emu_std_database
echo "Sequences path: $1"
# /mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass
echo "Output path: $2"
# /mnt/c/Users/Marcelo/Desktop/results
export EMU_DATABASE_DIR=$1
echo $EMU_DATABASE_DIR
# Get the directories for all the barcodes
export barcode_directories=$2/barcode*
echo $barcode_directories