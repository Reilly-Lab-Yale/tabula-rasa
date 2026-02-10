#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=day
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate env_tzinb

ipython fit.py

echo "EXITING SHELL"
