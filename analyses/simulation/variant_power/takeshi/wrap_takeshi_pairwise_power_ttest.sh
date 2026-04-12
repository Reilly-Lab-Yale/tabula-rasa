#!/bin/bash
#SBATCH --job-name=takeshi_pw
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module load miniconda
conda activate tz

python takeshi_pairwise_power_ttest.py all
echo "EXITING SHELL"
