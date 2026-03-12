#!/bin/bash
set -euo pipefail

#SBATCH --job-name=fit
#SBATCH --partition=ycga
#SBATCH --time=1-12:10:00
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module reset
module load miniconda
#conda activate env_tzinb #bouchet
conda activate env_tensorzinb #mccleary

ipython shendure_consider_missing.py

echo "EXITING SHELL"
