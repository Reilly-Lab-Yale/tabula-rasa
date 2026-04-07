#!/bin/bash
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH --job-name=seelig_bounds
#SBATCH -o seelig_bounds_%j.out

module load miniconda
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/model_fitting/bounds

python extract_seelig_bounds.py
