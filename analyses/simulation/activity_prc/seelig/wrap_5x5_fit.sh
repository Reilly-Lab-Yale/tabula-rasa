#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --gpus=h200:1
#SBATCH -c 8
#SBATCH --mem=256G
#SBATCH -t 2-00:00:00
#SBATCH -J seelig_5x5_fit
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Seelig CM: consider_missing expands data to full Cartesian product.
# 6 LocalCluster workers (processes=True) to parallelize 5 reps per GT draw.
# Empirical CM fit peaked at ~49GB RAM per ortho; 256GB for headroom.
# ETA: ~4.5h (5 GT draws x ~55min each with 5-way rep parallelism)

module load miniconda
conda activate tz
module load CUDA cuDNN

export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH
cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc/seelig

python seelig_5x5_activity.py fit
