#!/bin/bash
#SBATCH --job-name=shend_pow
#SBATCH --partition=ycga
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate tz

jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 \
    shendure_power_mwu_all_cell_types.ipynb

echo "EXITING SHELL"
