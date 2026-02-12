#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=week
#SBATCH --time=2-00:10:00
#SBATCH -c 6
#SBATCH --mem=512G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module reset
module load miniconda
#conda activate env_tzinb #bouchet
conda activate env_tensorzinb #mccleary

ipython fit.py

echo "EXITING SHELL"
