#!/bin/bash
#SBATCH --job-name=shend_pw
#SBATCH --partition=ycga
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm-pw-%j.out
#SBATCH --error=slurm-pw-%j.err

module load miniconda
conda activate tz

jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 \
    shendure_pairwise_power_mwu.ipynb

exit_code=$?

if [ $exit_code -eq 0 ]; then
    NTFY_TITLE="Shendure pairwise power done" NTFY_TAGS="white_check_mark" \
        notify-job "shendure_pairwise_power_mwu completed successfully on $(hostname)"
else
    NTFY_TITLE="Shendure pairwise power FAILED" NTFY_TAGS="warning,x" \
        notify-job "shendure_pairwise_power_mwu exited with code $exit_code on $(hostname). Check slurm-pw-${SLURM_JOB_ID}.err"
fi

echo "EXITING SHELL"
