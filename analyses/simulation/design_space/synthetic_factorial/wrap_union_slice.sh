#!/bin/bash
#SBATCH --job-name=synfac
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

# Slice-grain array wrapper for the union sweep.
# Each task is one slice (out of n_slices=50). Slice i processes samples
# where idx % 50 == i, so ~100 samples per slice sequentially.
# Slice driver reuses a single SLURMCluster across all its samples, so
# worker sbatches happen once (~5 per slice) instead of per-sample.
#
# Use:  sbatch --array=0-49 wrap_union_slice.sh

module load miniconda
conda activate tz
python synthetic_factorial.py simulate union "$SLURM_ARRAY_TASK_ID" 50
echo "EXITING SHELL"
