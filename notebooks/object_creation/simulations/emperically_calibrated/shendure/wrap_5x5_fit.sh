#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --gpus=h200:1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -J shend_5x5_fit
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load miniconda
conda activate tz
module load CUDA cuDNN

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/simulations/emperically_calibrated/shendure

python shendure_5x5_activity.py fit
