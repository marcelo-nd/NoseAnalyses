*add bin folder to PATH
if [ -d "$HOME/marcelo/.local/bin" ] ; then   PATH="$PATH:$HOME//marcelo/.local/bin"; fi

export EMU_DATABASE_DIR=./emu_std_database

export barcode1files=/mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass/barcode01/fastaq/*.fastq

for f in $barcode1files;
do emu abundance $f --keep-files --keep-counts --keep-read-assignments --output-dir ./results1 ;
done

export barcode2files=/mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass/barcode02/fastaq/*.fastq

for f in $barcode2files;
do emu abundance $f --keep-files --keep-counts --keep-read-assignments --output-dir ./results2 ;
done

export barcode3files=/mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass/barcode03/fastaq/*.fastq

for f in $barcode3files;
do emu abundance $f --keep-files --keep-counts --keep-read-assignments --output-dir ./results3 ;
done

export barcode4files=/mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass/barcode04/fastaq/*.fastq

for f in $barcode4files;
do emu abundance $f --keep-files --keep-counts --keep-read-assignments --output-dir ./results4 ;
done


# Create custom database
emu build-database LaCa --sequences /mnt/c/Users/Marcelo/Desktop/16s_nasal_strains.fasta --seq2tax /mnt/c/Users/Marcelo/Desktop/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/Marcelo/Desktop/taxonomy.tsv


 . /mnt/c/Users/Marcelo/Documents/GitHub/NoseGenomicAnalyses/EmuWrapper/emu_wrapper.sh ./emu_std_database /mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass/


 #fix line endings
  sed -i 's/\r//' /mnt/c/Users/Marcelo/Documents/GitHub/NoseGenomicAnalyses/EmuWrapper/emu_wrapper.sh