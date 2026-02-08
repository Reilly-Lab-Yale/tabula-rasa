#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=week
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=256G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate env_tzinb

ipython fit.py

echo "EXITING SHELL"