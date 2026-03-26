#!/bin/bash
#SBATCH --job-name=shend_cal
#SBATCH --partition=ycga
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm-cal-%j.out
#SBATCH --error=slurm-cal-%j.err

module load miniconda
conda activate tz

jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 \
    shendure_calibration_mwu_all_cell_types.ipynb

exit_code=$?

if [ $exit_code -eq 0 ]; then
    NTFY_TITLE="Shendure calibration done" NTFY_TAGS="white_check_mark" \
        notify-job "shendure_calibration_mwu completed successfully on $(hostname)"
else
    NTFY_TITLE="Shendure calibration FAILED" NTFY_TAGS="warning,x" \
        notify-job "shendure_calibration_mwu exited with code $exit_code on $(hostname). Check slurm-cal-${SLURM_JOB_ID}.err"
fi

echo "EXITING SHELL"
