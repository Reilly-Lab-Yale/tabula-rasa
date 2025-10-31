#!/bin/bash
#SBATCH -p day
#SBATCH -t 8:00:00
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -o %j_cpu2_dense
module load miniconda
conda activate env_tensorzinb_redo

python speedtest.py cpu2dense predense
