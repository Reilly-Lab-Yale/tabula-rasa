#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --job-name=phantom_obs_zinb_bct
#SBATCH --cpus-per-task=2
#SBATCH --output=job_bct_%j.out
#SBATCH --error=job_bct_%j.err

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

python3 fit_phantom_by_cell_type.py
