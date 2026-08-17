#!/bin/bash
#SBATCH --job-name=seelig_pow
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --exclude=a1132u18n02
module load miniconda
conda activate tz

python seelig_power_ttest_all_cell_types.py all
echo "EXITING SHELL"
