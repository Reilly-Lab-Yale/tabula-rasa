#!/bin/bash
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --job-name=seelig_collision
#SBATCH -o seelig_collision_%j.out

module load miniconda/24.11.3
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa

python analyses/model_fitting/bounds/seelig_collision_rate.py
