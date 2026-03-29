#!/bin/bash
#SBATCH -J shend_cm_nb_ct
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -t 2-23:40:00
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module reset
module load miniconda
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos/shendure_cm_nb
ipython fit_by_cell_type.py

echo "EXITING SHELL"
