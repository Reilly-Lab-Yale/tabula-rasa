#!/bin/bash
#SBATCH --job-name=prc_curves
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -uo pipefail
module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc
python plot_curves.py
echo "exit=$?"
