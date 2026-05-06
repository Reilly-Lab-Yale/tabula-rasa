#!/bin/bash
# launch.sh -- orchestrate a full synthetic-factorial run via N parallel slices.
#
# Steps:
#   1. Materialize the LHS samples table (synchronous, fast).
#   2. Submit N independent simulate-slice sbatch jobs, one per slice.
#      Each job has its own dask SLURMCluster, so there is no shared scheduler
#      and no nested-dask contention.
#   3. Submit a single plot-phase sbatch job with --dependency=afterok on all
#      slice jobs. The plot job aggregates outputs across slices.
#
# Usage:  ./launch.sh [mode]          mode: smoke | pilot | full   (default: full)
#
# Output: prints the slice job IDs and the plot job ID. Monitor with
#   squeue -u $USER

mode=${1:-full}

cd "$(dirname "$0")"

module load miniconda
conda activate tz

# Step 1: write the samples table so all slices share the same draw.
python synthetic_factorial.py samples "$mode"

# Step 2: read n_slices from the mode config and submit slice jobs.
n_slices=$(python -c "
import sys
sys.path.insert(0,'.')
import synthetic_factorial as s
print(s.MODES['$mode']['n_slices'])
")
echo "Submitting $n_slices simulate-slice jobs for mode=$mode"

slice_jids=""
for i in $(seq 0 $((n_slices-1))); do
    out=$(sbatch wrap_synthetic_factorial.sh simulate "$mode" "$i" "$n_slices")
    jid=$(echo "$out" | awk '{print $4}')
    slice_jids="${slice_jids}:${jid}"
    echo "  slice $i  jobid=$jid"
done
slice_jids=${slice_jids#:}  # strip leading colon

# Step 3: dependent plot job.
plot_out=$(sbatch --dependency=afterok:"$slice_jids" \
    wrap_synthetic_factorial.sh plot "$mode")
plot_jid=$(echo "$plot_out" | awk '{print $4}')
echo "Plot job $plot_jid will run after slices $slice_jids"
echo
echo "Monitor:  squeue -u \$USER"
echo "Slice IDs: $slice_jids"
echo "Plot ID:   $plot_jid"
