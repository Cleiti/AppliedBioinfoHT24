#!/bin/bash -l
#SBATCH -A naiss2023-22-1312
#SBATCH -p core
#SBATCH -n 4
#SBATCH -e error.log
#SBATCH -t 01:00:00
#SBATCH -J gape_blast
#SBATCH --mail-type=ALL
#SBATCH --mail-user gabrielandre.pettersson.9739@student.uu.se

#Load stuff

module load bioinfo-tools
module load blast/2.15.0+

# Create BLAST database for clavata
makeblastdb -in ../input/clavipes_protein.fasta -dbtype prot -out dbs/clavipes_db -parse_seqids

# Compare HP126 with reference genome using 99% identity threshold
blastp -query ../input/W1_S2_3830_longest_ORFs_aa.fasta -db clavipes_db -out ../output/clavipes_3830.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4
       
# Cut out the sequence ID for the matches

cut -f 2 ../output/clavipes_3830.txt | sort | uniq > ../input/matched_ref_clavipes.txt

# Tie the sequence ID to the database

blastdbcmd -db dbs/clavipes_db -entry_batch ../input/matched_ref_clavipes.txt -outfmt "%f" -out ../input/matched_ref_clavipes.fasta

# Reverse blast db

makeblastdb -in ../input/W1_S2_3830_longest_ORFs_aa.fasta -dbtype prot -out dbs/rev_clavipes_db

# Reverse blast

blastp -query ../input/matched_ref_clavipes.fasta -db dbs/rev_clavipes_db -out ../output/reverse_clavipes_3830.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4
       
awk 'NR==FNR {a[$1]=$2; next} ($1 in a) && a[$1]==$2' ../output/clavipes_3830.txt ../output/reverse_clavipes_3830.txt > ../output/reciprocal_hits_3830.txt