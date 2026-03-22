#!/bin/bash
#SBATCH --job-name=cohen_pow
#SBATCH -A prio_skr2
#SBATCH --partition=priority
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module load miniconda
conda activate env_tzinb

jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 \
    cohen_power_mwu_all_cell_types.ipynb
echo "EXITING SHELL"
