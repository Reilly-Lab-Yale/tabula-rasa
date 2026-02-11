#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=ycga_long
#SBATCH --time=5-00:10:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=239G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module reset
module load miniconda
#conda activate env_tzinb #bouchet
conda activate env_tensorzinb #mccleary

ipython fit.py

echo "EXITING SHELL"
