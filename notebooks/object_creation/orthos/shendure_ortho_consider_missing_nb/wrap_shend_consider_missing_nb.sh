#!/bin/bash
#SBATCH --job-name=fit
#SBATCH --partition=ycga
#SBATCH --time=2-23:50:00
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module reset
module load miniconda
#conda activate env_tzinb #bouchet
#conda activate env_tensorzinb #mccleary
conda activate tz #mccleary (sparse TensorZINB)

ipython shendure_consider_missing_nb.py

echo "EXITING SHELL"
