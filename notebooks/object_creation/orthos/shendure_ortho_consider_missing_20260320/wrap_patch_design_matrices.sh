#!/bin/bash
#SBATCH --job-name=patch_design
#SBATCH --partition=ycga
#SBATCH --time=2-23:50:00
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module reset
module load miniconda
conda activate tz

ipython patch_design_matrices.py

echo "EXITING SHELL"
