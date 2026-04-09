#!/bin/bash
#SBATCH --job-name=shendure_pow
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module load miniconda
conda activate tz

python shendure_power_ttest_all_cell_types.py all
echo "EXITING SHELL"
