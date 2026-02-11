#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=ycga
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module reset
module load miniconda
#conda activate env_tzinb #bouchet
conda activate env_tensorzinb #mccleary

ipython fit.py

echo "EXITING SHELL"
