#!/bin/bash
#SBATCH --job-name=seelig_replot
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -uo pipefail
module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/variant_power/seelig
python seelig_pairwise_power_mwu.py plot
echo "exit=$?"
