#!/bin/bash
#SBATCH -J shend_cm_nb_cre
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -t 2-23:50:00
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module reset
module load miniconda
conda activate tz
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos/shendure_cm_nb
ipython fit_by_cre.py

echo "EXITING SHELL"
