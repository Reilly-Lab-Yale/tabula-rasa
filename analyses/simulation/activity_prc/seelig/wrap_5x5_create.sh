#!/bin/bash
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -J seelig_5x5_create
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load miniconda
conda activate tz

export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc/seelig

python seelig_5x5_activity.py create
