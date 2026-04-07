#!/bin/bash
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -c 2
#SBATCH --mem=128G
#SBATCH -t 4:00:00
#SBATCH -J cohen_5x5_create
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load miniconda
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc/cohen

python cohen_5x5_activity.py create
