#!/bin/bash
# Resubmit the 6 slices that died on 2026-05-05 from sbatch rate-limit storm.
# Slices 0,1,2,3,4,9 had zero-worker dask clusters and failed every sample.
# Slices 5,6,7,8 (jobs 10768460-10768463) are still running healthily.
# Plot job 10768465 is pending; we update its --dependency afterward.
#
# Stagger sbatches by 60s each to avoid re-tripping the rate limit and to
# let cluster.adapt() in each new driver kick in gradually.

set -u
mode=full
plot_jid=10768465
healthy="10768460:10768461:10768462:10768463"

cd "$(dirname "$0")"

new_jids=""
for slice in 0 1 2 3 4 9; do
    out=$(sbatch wrap_synthetic_factorial.sh simulate "$mode" "$slice" 10)
    jid=$(echo "$out" | awk '{print $4}')
    echo "  slice $slice  jobid=$jid"
    new_jids="${new_jids}:${jid}"
    sleep 60
done
new_jids=${new_jids#:}

# Update plot job dependency to include the new slice IDs.
all_deps="${healthy}:${new_jids}"
echo "Updating plot job $plot_jid dependency to: afterok:$all_deps"
scontrol update jobid=$plot_jid Dependency=afterok:$all_deps

echo
echo "New slice IDs: $new_jids"
echo "Plot $plot_jid will fire after: $all_deps"
