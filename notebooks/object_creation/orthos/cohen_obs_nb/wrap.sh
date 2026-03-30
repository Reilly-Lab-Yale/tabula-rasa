#!/bin/bash
#SBATCH -J cohen_obs_nb
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module reset
module load miniconda
conda activate tz
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos/cohen_obs_nb
ipython fit.py

echo "EXITING SHELL"
