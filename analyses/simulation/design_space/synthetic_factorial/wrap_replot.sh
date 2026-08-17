#!/bin/bash
#SBATCH --job-name=ds_replot
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -uo pipefail
module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/design_space/synthetic_factorial
python synthetic_factorial.py replot union
echo "exit=$?"
