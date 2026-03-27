#!/bin/bash
#SBATCH -J seelig_ct
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -t 2-00:00:00
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module reset
module load miniconda
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos
ipython fit_seelig_ortho_by_cell_type.py

echo "EXITING SHELL"
