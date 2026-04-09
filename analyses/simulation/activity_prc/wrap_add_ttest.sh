#!/bin/bash
#SBATCH --job-name=add_counts_tests
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=priority
#SBATCH --account=prio_skr2

module load miniconda
conda activate tz

export PYTHONPATH=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa:$PYTHONPATH

python add_counts_tests.py all
