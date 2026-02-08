#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=week
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate env_tzinb

ipython fit.py

echo "EXITING SHELL"