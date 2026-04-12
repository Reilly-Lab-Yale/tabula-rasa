#!/bin/bash
#SBATCH --job-name=shendure_mwucal
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --exclude=a1132u18n02,a1130u22n02
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
module load miniconda
conda activate tz

python shendure_calibration_mwu_followup.py all
echo "EXITING SHELL"
