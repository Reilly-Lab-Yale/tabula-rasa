#!/bin/bash
#SBATCH --job-name=fit
#SBATCH -A prio_reilly
#SBATCH --partition=priority
#SBATCH --time=6-00:10:00
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
