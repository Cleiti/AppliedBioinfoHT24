#!/bin/bash -l
#SBATCH -A naiss2023-22-1312
#SBATCH -p core
#SBATCH -n 4
#SBATCH -e error.log
#SBATCH -t 06:00:00
#SBATCH -J gape_blast
#SBATCH --mail-type=ALL
#SBATCH --mail-user gabrielandre.pettersson.9739@student.uu.se

#Load stuff

module load bioinfo-tools
module load blast/2.15.0+

# Create BLAST database for clavata
makeblastdb -in ../input/clavata_protein.fasta -dbtype prot -out clavata_db



# Compare HP126 with reference genome using 99% identity threshold
blastp -query ../input/W1_S5_4676_longest_ORFs_aa.fasta -db clavata_db -out ../output/clavata_4676.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4

# Compare DV3 with reference genome using 99% identity threshold
blastp -query ../input/W1_S8_4682_longest_ORFs_aa.fasta -db clavata_db -out ../output/clavata_4682.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4

# Compare genome1 with genome2 using 99% identity threshold
blastp -query ../input/W1_S9_4683_longest_ORFs_aa.fasta -db clavata_db -out ../output/clavata_4683.txt \
       -outfmt 6 -evalue 1e-5 -num_threads 4