#!/bin/bash
#SBATCH --job-name=seelig_pw
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --exclude=a1132u18n02,a1130u22n02
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate tz

python seelig_pairwise_power_mwu.py all
exit_code=$?

if [ $exit_code -eq 0 ]; then
    NTFY_TITLE="Seelig pairwise power done" NTFY_TAGS="white_check_mark" \
        notify-job "seelig_pairwise_power_mwu completed successfully on $(hostname)"
else
    NTFY_TITLE="Seelig pairwise power FAILED" NTFY_TAGS="warning,x" \
        notify-job "seelig_pairwise_power_mwu exited with code $exit_code on $(hostname). Check slurm-${SLURM_JOB_ID}.err"
fi

echo "EXITING SHELL"
