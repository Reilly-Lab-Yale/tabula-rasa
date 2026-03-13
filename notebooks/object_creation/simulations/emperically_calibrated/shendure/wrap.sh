#!/bin/bash
#SBATCH --job-name=shend_pow
#SBATCH -A prio_reilly
#SBATCH --partition=priority
#SBATCH --time=18:10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate env_tzinb

jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 \
    shendure_power_mwu_all_cell_types.ipynb

echo "EXITING SHELL"
