#!/bin/bash
#SBATCH --job-name=dsweep
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Usage: sbatch wrap_design_space_sweep.sh [phase] [axis] [test_type]
#   phase:      simulate | plot | all  (default: all)
#   axis:       bcs_per_cre | n_cells | moi | n_cres  (default: bcs_per_cre)
#   test_type:  ttest | mwu  (default: ttest, only used in plot phase)

module load miniconda
conda activate tz

python design_space_sweep.py "$@"
echo "EXITING SHELL"
