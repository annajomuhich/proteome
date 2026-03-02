#!/bin/bash
#SBATCH -D /home/ajmuhich/
#SBATCH -o /home/ajmuhich/slurm-log/tc_at_model_timepoint_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/tc_at_model_timepoint_stderr-%j.txt
#SBATCH -J tc_at_model_timepoint
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
COUNTS_FILE1=$1
FAILED_GENES_FILE=$2    
OUTPUT_DIR=$3

Rscript ~/proteome/tc_at_model_timepoint.R \
"$COUNTS_FILE1" \
"$FAILED_GENES_FILE" \
"$OUTPUT_DIR"
