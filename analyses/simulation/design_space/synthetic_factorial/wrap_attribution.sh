#!/bin/bash
#SBATCH --job-name=synfac_attr
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Regenerate the 3-pair per-axis attribution figures from the union LHS
# (output/samples_power_union.parquet). Single process, 3 pairs run
# sequentially; LOESS bootstrap is single-threaded so a few CPUs suffice.
# Off the login node to avoid hammering the shared 1-CPU cgroup.

module load miniconda
conda activate tz

# Do NOT cd "$(dirname "$0")" -- in slurm $0 is the staged copy in
# /var/spool/slurmd. Slurm sets WorkDir to the sbatch invocation cwd.

python attribution.py
echo "EXITING SHELL"
