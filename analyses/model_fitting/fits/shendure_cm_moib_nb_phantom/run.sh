#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=256G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=shend_cm_moib
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --exclude=a1132u18n02

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

python3 fit.py
