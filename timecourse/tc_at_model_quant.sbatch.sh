#!/bin/bash
#SBATCH -D /home/ajmuhich/
#SBATCH -o /home/ajmuhich/slurm-log/tc_at_model_quant_stdout-%j.txt
#SBATCH -e /home/ajmuhich/slurm-log/tc_at_model_quant_stderr-%j.txt
#SBATCH -J tc_at_model_quant
#SBATCH -t 2-00:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ajmuhich@ucdavis.edu

module load R

# Define input file paths
COUNTS_FILE1=$1
OUTPUT_DIR=$2

Rscript ~/proteome/tc_at_model_quant.R \
"$COUNTS_FILE1" \
"$OUTPUT_DIR"
