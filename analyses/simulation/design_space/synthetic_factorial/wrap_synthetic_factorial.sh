#!/bin/bash
#SBATCH --job-name=synfac
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Usage:
#   sbatch wrap_synthetic_factorial.sh <phase> <mode> [slice n_slices]
#     phase:    samples | simulate | plot | all
#     mode:     smoke | pilot | full
#     slice:    0..n_slices-1 (default 0; only meaningful for simulate)
#     n_slices: total slice count (default from MODES[mode]["n_slices"])
#
# Recommended: use ./launch.sh <mode> instead, which orchestrates samples
# materialization, N parallel simulate-slice jobs, and a dependent plot job.
#
# To re-plot from a relocated sim tree (e.g. cold storage):
#   SYNTHETIC_FACTORIAL_SIM_ROOT=/path/to/cold/synthetic_factorial_sims/<dated_dir> \
#     sbatch wrap_synthetic_factorial.sh plot full

module load miniconda
conda activate tz

python synthetic_factorial.py "$@"
echo "EXITING SHELL"
