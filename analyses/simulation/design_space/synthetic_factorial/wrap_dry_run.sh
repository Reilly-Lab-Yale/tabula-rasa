#!/bin/bash
#SBATCH --job-name=dryrun
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate tz
python dry_run_cohen_scale.py
echo "EXITING SHELL"
