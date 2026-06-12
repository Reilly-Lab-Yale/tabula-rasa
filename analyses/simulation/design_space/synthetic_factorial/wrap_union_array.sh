#!/bin/bash
#SBATCH --job-name=synfac
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

# Array wrapper for the union sweep. Each task processes one sample
# (slice_idx = SLURM_ARRAY_TASK_ID, n_slices = 5000). Use this with
#   sbatch --array=<lo>-<hi>%<concurrency> wrap_union_array.sh
# to avoid the per-hour sbatch rate limit that one-by-one submission
# tripped on 2026-05-11.

module load miniconda
conda activate tz
# Do NOT cd "$(dirname "$0")" -- in array mode $0 is slurm's staged copy
# in /var/spool, not the project location. Slurm sets WorkDir to the sbatch
# invocation cwd, which is what we want.
python synthetic_factorial.py simulate union "$SLURM_ARRAY_TASK_ID" 5000
echo "EXITING SHELL"
