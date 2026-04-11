#!/bin/bash
#SBATCH --job-name=seelig_cal
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --exclude=a1132u18n02
module load miniconda
conda activate tz

python seelig_calibration_ttest_all_cell_types.py all
echo "EXITING SHELL"
