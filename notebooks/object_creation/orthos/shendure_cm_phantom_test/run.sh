#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=phantom_cm_test
#SBATCH --cpus-per-task=2
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

python3 test_by_cre.py
