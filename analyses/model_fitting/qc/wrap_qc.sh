#!/bin/bash
#SBATCH -J ortho_qc
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH -t 4:00:00
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -o slurm-qc-%j.out
#SBATCH -e slurm-qc-%j.err

# Usage: sbatch wrap_qc.sh <dataset> <ortho_name>
# Example: sbatch wrap_qc.sh shendure shendure_obs_nb_20260329

module reset
module load miniconda
conda activate tz
export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/model_fitting/qc
ipython run_qc.py -- --dataset "$1" --ortho "$2"

echo "EXITING SHELL"
