#!/bin/bash
#SBATCH -J sync_pow_sims
#SBATCH -p day
#SBATCH -t 12:00:00
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -o slurm-sync-%j.out
#SBATCH -e slurm-sync-%j.err

set -euo pipefail

SRC="bouchet.ycrc.yale.edu:/nfs/roberts/project/pi_skr2/shared/tabula_data"
DST="/gpfs/gibbs/pi/reilly/tabula_data"

echo "Starting sync at $(date)"

rsync -av --progress \
    "${SRC}/2026-03-18_shendure_pow" \
    "${DST}/"

rsync -av --progress \
    "${SRC}/2026-03-21_cohen_pow" \
    "${DST}/"

echo "Sync complete at $(date)"
