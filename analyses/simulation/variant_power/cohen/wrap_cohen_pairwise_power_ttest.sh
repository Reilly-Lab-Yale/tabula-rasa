#!/bin/bash
#SBATCH --job-name=cohen_pw_tt
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --exclude=a1132u18n02
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
conda activate tz

python cohen_pairwise_power_ttest.py all
exit_code=$?

if [ $exit_code -eq 0 ]; then
    NTFY_TITLE="Cohen pairwise power (t-test) done" NTFY_TAGS="white_check_mark" \
        notify-job "cohen_pairwise_power_ttest completed successfully on $(hostname)"
else
    NTFY_TITLE="Cohen pairwise power (t-test) FAILED" NTFY_TAGS="warning,x" \
        notify-job "cohen_pairwise_power_ttest exited with code $exit_code on $(hostname). Check slurm-${SLURM_JOB_ID}.err"
fi

echo "EXITING SHELL"
