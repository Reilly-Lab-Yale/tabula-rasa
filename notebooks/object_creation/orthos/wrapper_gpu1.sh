#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH -t 8:00:00
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -o %j_gpu1_sparse
module load miniconda
conda activate env_tensorzinb_redo

python speedtest.py gpu1sparse sparse

#SBATCH --partition=gpu_h200
