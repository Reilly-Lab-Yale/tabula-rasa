#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=day
#SBATCH --time=08:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate env_tzinb

ipython fit.py

echo "EXITING SHELL"