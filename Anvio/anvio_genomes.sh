export user_path=/mnt/e/KarosGenomes

cd $user_path

mkdir $user_path/anvio_analyses/contigs-fasta_files

# Reformat fasta files to contigs-fasta files (anvio-style fasta).
for file in *.{fna,fasta}
do
    printf "%s \n" "$file"
    anvi-script-reformat-fasta $file\
    -o "$user_path/anvio_analyses/contigs-fasta_files/${file}"\
    --simplify-names
done

cd $user_path/anvio_analyses/contigs-fasta_files

mkdir $user_path/anvio_analyses/contigs-db

# Generate contigs-db files (entry-point files for anvio).
for file in *.fna
do
    printf "%s \n" "$file"
    anvi-gen-contigs-database -f $file\
    -o "$user_path/anvio_analyses/contigs-db/${file/%fna/db}"\

done

cd $user_path/anvio_analyses/contigs-db

# Run HMMs (annotate your genes against an hmm-source).
for file in *.db
do
    printf "%s \n" "$file"
    anvi-run-hmms -c $file\
    --num-threads 6
done

anvi-setup-scg-taxonomy

# Run SCG Taxonomy, associates its single-copy core gene with taxonomic data.
for file in *.db
do
    printf "%s \n" "$file"
    anvi-run-scg-taxonomy -c $file\
    --min-percent-identity 90
done

# Scan tRNAs. Identify and store tRNA genes in a contigs database.
for file in *.db
do
    printf "%s \n" "$file"
    anvi-scan-trnas -c $file\

done

# Set up COG data
anvi-setup-ncbi-cogs

# Functional annotation with the COGs database.
for file in *.db
do
    printf "%s \n" "$file"
    anvi-run-ncbi-cogs -c $file\

done

# Set up KEGG data
anvi-setup-kegg-data

# Functional annotation with KEGG KOfam database.
for file in *.db
do
    printf "%s \n" "$file"
    anvi-run-kegg-kofams -c $file\
    -T 6

done

# Generate genomes-storage-db file. This is the input for the pangenome analyses.
anvi-gen-genomes-storage -e external_genomes.txt -o genomes-storage-db-GENOMES.db

#NOSE1
# Run Pangenome Analysis.
anvi-pan-genome -g genomes-storage-db-GENOMES.db \
--project-name "Nose_Pan" \
--output-dir NOSE \
--num-threads 6

# Displaying the pan genome
anvi-display-pan -p ./NOSE/Nose_Pan-PAN.db -g genomes-storage-db-GENOMES.db

#NOSE2
# Run Pangenome Analysis.
anvi-pan-genome -g genomes-storage-db-GENOMES.db \
--project-name "Nose_Pan" \
--output-dir NOSE2 \
--enforce-hierarchical-clustering \
--num-threads 6

# Displaying the pan genome
anvi-display-pan -p ./NOSE2/Nose_Pan-PAN.db -g genomes-storage-db-GENOMES.db

# NOSE3
# Run Pangenome Analysis.
anvi-pan-genome -g genomes-storage-db-GENOMES.db \
--project-name "Nose_Pan" \
--output-dir NOSE3 \
--min-occurrence 2 \
--num-threads 6

# Displaying the pan genome
anvi-display-pan -p ./NOSE3/Nose_Pan-PAN.db -g genomes-storage-db-GENOMES.db
