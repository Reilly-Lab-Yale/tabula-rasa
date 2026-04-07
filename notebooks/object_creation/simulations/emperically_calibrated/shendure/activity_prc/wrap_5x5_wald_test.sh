#!/bin/bash
#SBATCH -p priority_gpu
#SBATCH -A prio_skr2
#SBATCH --gpus=h200:1
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 4:00:00
#SBATCH -J shend_5x5_wald
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load miniconda
conda activate tz
module load CUDA cuDNN

cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/simulations/emperically_calibrated/shendure

python shendure_5x5_activity.py wald_precomp
python shendure_5x5_activity.py test
python shendure_5x5_activity.py metrics
