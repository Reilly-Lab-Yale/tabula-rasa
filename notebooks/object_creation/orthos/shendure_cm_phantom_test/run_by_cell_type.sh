#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=phantom_cm_bct
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

python3 test_by_cell_type.py
