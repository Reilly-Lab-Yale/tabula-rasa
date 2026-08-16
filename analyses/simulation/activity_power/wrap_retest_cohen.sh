#!/bin/bash
#SBATCH --job-name=cohen_retest
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -uo pipefail
module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_power
python retest_arms.py \
  /nfs/roberts/scratch/pi_skr2/mcn26/restored/2026-04-08_cohen_pow \
  cohen_power_df_2026-08-13 \
  --arms mwu,mwu_deflated
echo "exit=$?"
