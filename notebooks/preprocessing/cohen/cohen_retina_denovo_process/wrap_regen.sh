#!/bin/bash
#SBATCH -J cohen_regen
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_regen_%j.out
#SBATCH --error=slurm_regen_%j.err

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/preprocessing/cohen/cohen_retina_denovo_process

module reset
module load miniconda
conda activate tz
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH

python regen_scmpra_object.py
