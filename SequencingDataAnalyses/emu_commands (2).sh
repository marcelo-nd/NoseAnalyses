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
emu build-database LaCa --sequences /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/LaCa_sequences.fasta --seq2tax /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/LaCa_taxonomy.tsv


 . /mnt/c/Users/Marcelo/Documents/GitHub/NoseGenomicAnalyses/EmuWrapper/emu_wrapper.sh -d ./emu_std_database -s /mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass/


 #fix line endings
  sed -i 's/\r//' /mnt/c/Users/Marcelo/Documents/GitHub/EmuWrapper/emu_wrapper.sh

# ./emu_std_database
# ./LaCa

# /mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass

# /mnt/c/Users/Marcelo/Desktop/results


. emu_wrapper.sh -d /home/marcelo/LaCa -s /mnt/c/Users/Marcelo/Desktop/230815_Marcelo/230816-Marcelo/230816-Marcelo/20230816_1015_MN31656_FAR91361_edfa8c58/fastq_pass -o /mnt/c/Users/Marcelo/Desktop/results -b "all" -n "test"


############

conda activate emu_py37
 
cd /mnt/c/Users/Marcelo\ Navarro/Documents/GitHub/EmuWrapper
 

###########
# Create custom database
emu build-database LaCa --sequences /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/LaCa_sequences.fasta --seq2tax /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/LaCa_seq2tax.map --taxonomy-list /mnt/c/Users/'Marcelo Navarro'/'OneDrive - UT Cloud'/'NoseSynCom Project'/LaCa_16s_emu_database/LaCa_taxonomy.tsv

conda install r

sudo apt-get install libhdf5-dev

cd /mnt/c/Users/'Marcelo Navarro'/Documents/GitHub/EmuWrapper

. emu_wrapper.sh -d /home/marcelo/LaCa -s /mnt/e/NoseSynComProject/SequencingData/16S_First_Step_17_01_24/no_sample/20240117_1510_MN45148_FAY00412_e8d22899/fastq_pass -o /mnt/e/NoseSynComProject/SequencingData/16S_First_Step_17_01_24/results -b "all" -n "Test"


. emu_wrapper.sh -d /home/marcelo/LaCa -s /mnt/e/NoseSynComProject/SequencingData/Test2/no_sample/20240124_1634_MN45148_FAY00412_07f41999/fastq_pass -o /mnt/e/NoseSynComProject/SequencingData//Test2/results -b "all" -n "Test"

. emu_wrapper.sh -d test -s test -o test -b "test" -n "test"

. emu_wrapper.sh -d /home/marcelo/LaCa -s /mnt/e/NoseSynComProject/SequencingData/Test2/no_sample/20240124_1634_MN45148_FAY00412_07f41999/fastq_pass -o /mnt/e/NoseSynComProject/SequencingData/Test2/results



. emu_wrapper.sh -d /home/marcelo/LaCa -s /mnt/e/NoseSynComProject/SequencingData/SynComTest3_2/no_sample/20240212_1656_MN45148_FAY00412_e1705dac/fastq_pass -o /mnt/e/NoseSynComProject/SequencingData/SynComTest3_2/results
