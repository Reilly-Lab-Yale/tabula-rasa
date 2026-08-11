#!/bin/bash
#SBATCH --job-name=activity_volcano
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module load miniconda
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/empirical_testing/activity_volcano
python activity_volcano.py
echo "EXITING SHELL"
