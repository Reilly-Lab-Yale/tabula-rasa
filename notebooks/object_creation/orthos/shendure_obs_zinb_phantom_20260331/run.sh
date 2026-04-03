#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=phantom_obs_zinb
#SBATCH --cpus-per-task=2
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

python3 fit_phantom.py
