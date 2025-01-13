#!/bin/bash -l
#SBATCH -A naiss2023-22-1312
#SBATCH -p core
#SBATCH -n 4
#SBATCH -e error.log
#SBATCH -t 02:00:00
#SBATCH -J gape_rec_blast
#SBATCH --mail-type=ALL
#SBATCH --mail-user gabrielandre.pettersson.9739@student.uu.se

#Load stuff

module load bioinfo-tools
module load blast/2.15.0+

# Create BLAST database for clavata
makeblastdb -in ../input/clavata_protein.fasta -dbtype prot -out dbs/clavata_db -parse_seqids

# Cut out the sequence ID for the matches

cut -f 2 ../output/clavata_4676.txt | sort | uniq > ../input/matched_ref_4676.txt

cut -f 2 ../output/clavata_4682.txt | sort | uniq > ../input/matched_ref_4682.txt

cut -f 2 ../output/clavata_4683.txt | sort | uniq > ../input/matched_ref_4683.txt


# Tie the sequence ID to the database

blastdbcmd -db dbs/clavata_db -entry_batch ../input/matched_ref_4676.txt -outfmt "%f" -out ../input/matched_ref_4676.fasta
blastdbcmd -db dbs/clavata_db -entry_batch ../input/matched_ref_4682.txt -outfmt "%f" -out ../input/matched_ref_4682.fasta
blastdbcmd -db dbs/clavata_db -entry_batch ../input/matched_ref_4683.txt -outfmt "%f" -out ../input/matched_ref_4683.fasta

# Reverse blast db

makeblastdb -in ../input/W1_S5_4676_longest_ORFs_aa.fasta -dbtype prot -out dbs/rev_4676_db

makeblastdb -in ../input/W1_S8_4682_longest_ORFs_aa.fasta -dbtype prot -out dbs/rev_4682_db

makeblastdb -in ../input/W1_S9_4683_longest_ORFs_aa.fasta -dbtype prot -out dbs/rev_4683_db

# Reverse blast

blastp -query ../input/matched_ref_4676.fasta -db dbs/rev_4676_db -out ../output/reverse_clavata_4676.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4


blastp -query ../input/matched_ref_4682.fasta -db dbs/rev_4682_db -out ../output/reverse_clavata_4682.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4


blastp -query ../input/matched_ref_4683.fasta -db dbs/rev_4683_db -out ../output/reverse_clavata_4683.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4
       
       
awk 'NR==FNR {a[$1]=$2; next} ($1 in a) && a[$1]==$2' ../output/clavata_4676.txt ../output/reverse_clavata_4676.txt > ../output/reciprocal_hits_4676.txt
       
awk 'NR==FNR {a[$1]=$2; next} ($1 in a) && a[$1]==$2' ../output/clavata_4682.txt ../output/reverse_clavata_4682.txt > ../output/reciprocal_hits_4682.txt
 
awk 'NR==FNR {a[$1]=$2; next} ($1 in a) && a[$1]==$2' ../output/clavata_4683.txt ../output/reverse_clavata_4683.txt > ../output/reciprocal_hits_4683.txt      
       