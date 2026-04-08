#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --gpus=h200:1
#SBATCH -c 8
#SBATCH --mem=256G
#SBATCH -t 3-00:00:00
#SBATCH -J seelig_5x5_wald
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Seelig Wald: 1344-column by_cell_type models on CM-expanded data.
# ~3000s/model on H200, 2 models per ortho, 5 reps parallel per GT draw.
# ETA: ~8h (5 GT draws x ~100min each with 5-way rep parallelism)

module load miniconda
conda activate tz
module load CUDA cuDNN

export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc/seelig

python seelig_5x5_activity.py wald_precomp
python seelig_5x5_activity.py test
python seelig_5x5_activity.py metrics
