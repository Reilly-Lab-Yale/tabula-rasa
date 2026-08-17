#!/bin/bash
#SBATCH --job-name=shend_replot
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -uo pipefail
module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/variant_power/shendure
python replot_panels_mwu.py
echo "exit=$?"
