#!/bin/bash
# Rebuild the extended-grid pairwise-power aggregates from the 2026-08-13 sims.
# Read-only over the sim trees: the "plot" phase walks existing hs_pairwise
# results and writes *_pairwise_power_df.parquet. It never simulates, and it
# uses a LocalCluster, so it spawns no SLURM workers of its own.
#SBATCH --job-name=pw_aggregate
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -uo pipefail
module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

B=/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/analyses/simulation/variant_power
rc_all=0
for d in cohen shendure seelig; do
    echo "=============== ${d} plot phase ==============="
    cd "$B/$d" || { echo "no dir for $d"; rc_all=1; continue; }
    python "${d}_pairwise_power_mwu.py" plot
    rc=$?
    echo "${d} exit=${rc}"
    [ $rc -ne 0 ] && rc_all=$rc
done
echo "ALL DONE rc=${rc_all}"
exit $rc_all
