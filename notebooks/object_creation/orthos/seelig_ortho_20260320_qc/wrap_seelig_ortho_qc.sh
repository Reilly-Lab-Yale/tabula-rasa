#!/bin/bash
#SBATCH -p priority
#SBATCH -A prio_skr2
#SBATCH --job-name=seelig_ortho_qc
#SBATCH -t 02:00:00
#SBATCH -c 4
#SBATCH --mem=220G
#SBATCH --output=%x-%j.out

QC_DIR=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos/seelig_ortho_20260320_qc
cd "$QC_DIR"
/home/mcn26/.conda/envs/tz/bin/ipython "$QC_DIR/run_qc.py"
