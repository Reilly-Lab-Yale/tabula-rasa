#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --gpus=h200:1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -J seelig_5x5
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load miniconda
conda activate tz
module load CUDA cuDNN

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/activity_prc/seelig

python seelig_5x5_activity.py create
python seelig_5x5_activity.py fit
