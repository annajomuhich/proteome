#!/bin/bash
#SBATCH -D /home/ajmuhich/
#SBATCH -o /home/ajmuhich/slurm-log/mut_bc_model_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/mut_bc_model_stderr-%j.txt
#SBATCH -J mut_bc_model
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
COUNTS_FILE1=$1
OUTPUT_DIR=$2

Rscript ~/proteome/mut_bc_model.R \
"$COUNTS_FILE1" \
"$OUTPUT_DIR"
