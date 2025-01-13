#SBATCH -A naiss2023-22-1312          # Project allocation
#SBATCH -p core                       # Partition/queue
#SBATCH -n 1                          # Number of cores per job
#SBATCH -t 36:00:00                   # Maximum runtime
#SBATCH -J busco_species              # Job name
#SBATCH -e busco_species-%A_%a.err    # Error log for each array job
#SBATCH -o busco_species-%A_%a.out    # Output log for each array job
#SBATCH --array=0-225     
#SBATCH --mail-type=All
#SBATCH --mail-user=victor.englof.5352@student.uu.se


module load bioinfo-tools
module load BUSCO/4.1.4


source $AUGUSTUS_CONFIG_COPY

#directory containing the FASTA files
FASTA_DIR="/proj/uppstore2019013/nobackup/private/1000spider_master_project/fastafiles_renamed/renamed"

#List all FASTA files and store them in an array
FASTA_FILES=($(ls "$FASTA_DIR"/*.fasta))

#Get the file corresponding to this job's array index
FASTA_FILE=${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}

#Extract the species name (basename without .fasta extension)
SPECIES_NAME=$(basename "$FASTA_FILE" .fasta)

#Define the output directory based on species name
OUTPUT_DIR="/proj/uppstore2019013/nobackup/private/1000spider_master_project/busco/busco_results/${SPECIES_NAME}_busco_res"

#Create output directory if it doesn't already exist
mkdir -p "$OUTPUT_DIR"

run_BUSCO.py -i "$FASTA_FILE" -o "$OUTPUT_DIR" -l $BUSCO_LINEAGE_SETS/arachnida_odb10 -m prot