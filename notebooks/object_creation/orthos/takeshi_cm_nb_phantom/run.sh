#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --job-name=tak_cm_nb
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

python3 fit.py
