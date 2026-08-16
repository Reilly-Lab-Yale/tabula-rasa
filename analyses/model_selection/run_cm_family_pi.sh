#!/bin/bash
#SBATCH --job-name=cm_family_pi
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module load miniconda
conda activate tz

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/model_selection
python cm_family_pi.py
echo "EXITING SHELL"
