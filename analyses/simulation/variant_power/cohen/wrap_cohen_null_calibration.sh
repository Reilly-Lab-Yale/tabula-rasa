#!/bin/bash
#SBATCH --job-name=cohen_null
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --exclude=a1132u18n02
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate tz

python cohen_pairwise_null_calibration_ttest.py all
exit_code=$?

if [ $exit_code -eq 0 ]; then
    NTFY_TITLE="Cohen null calibration done" NTFY_TAGS="white_check_mark" \
        notify-job "cohen_pairwise_null_calibration_ttest completed successfully on $(hostname)"
else
    NTFY_TITLE="Cohen null calibration FAILED" NTFY_TAGS="warning,x" \
        notify-job "cohen_pairwise_null_calibration_ttest exited with code $exit_code on $(hostname). Check slurm-${SLURM_JOB_ID}.err"
fi

echo "EXITING SHELL"
