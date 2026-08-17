#!/bin/bash
#SBATCH --job-name=plot_curves
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module load miniconda
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc
python plot_curves.py
echo "EXITING SHELL"
