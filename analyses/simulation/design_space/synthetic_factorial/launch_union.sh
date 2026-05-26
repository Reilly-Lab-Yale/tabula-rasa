#!/bin/bash
# launch_union.sh -- submit sample-grain drivers for the union sweep.
#
# 5000 single-sample drivers (one per LHS draw). Each driver runs all 5
# library reps of its one sample, then exits. Inline cache+prune in
# _run_one_sample keeps peak scratch usage bounded (~10 GB/sample x
# concurrent sample count, then collapsed to ~10 KB).
#
# Sample-grain instead of slice-grain because the heaviest samples need
# >24h of driver wall time when run sequentially. One sample per driver
# fits comfortably in the 24h limit.
#
# Submits in waves to avoid sbatch rate limits (the original 2026-05-05
# pilot tripped the YCRC per-hour limit with a single burst).
#
# Usage:  ./launch_union.sh [wave_size] [wave_delay_seconds]
#   wave_size:        sbatches per wave (default 100)
#   wave_delay_seconds: pause between waves (default 60)
#
# Output: jid log printed line-per-sample; full list saved to
#   launch_union_<timestamp>.jids for later monitoring.

# Note: NOT using `set -e` because conda activate's internal helpers
# emit non-zero exits during normal operation, which would abort the
# launcher before the first sbatch.
set -uo pipefail

wave_size=${1:-100}
wave_delay=${2:-60}

cd "$(dirname "$0")"
module load miniconda
conda activate tz

mode=union

# Step 1: write samples table (idempotent).
python synthetic_factorial.py samples "$mode"

# Step 2: enumerate sample row indices that still need work.
# A sample is "done" if its dir contains cached_results.parquet (post-prune)
# or N_LIBRARY_REPS valid completed reps (pre-prune).
to_do_file=$(mktemp)
python3 - << 'PY' > "$to_do_file"
import sys
sys.path.insert(0, ".")
import synthetic_factorial as sf
sf._apply_mode("union"); sf._apply_mode_paths("union")
samples = sf.get_samples()
n_total = len(samples)
n_lhs = int(sf.N_LHS)
n_reps = int(sf.N_LIBRARY_REPS)
todo = []
for i, row in samples.iterrows():
    sid = row["sample_id"]
    pt = sf._sample_dir(sid)
    if (pt / "cached_results.parquet").exists():
        continue
    if pt.exists() and sf._count_done(pt) >= n_reps:
        continue
    todo.append(i)
print(f"# total={n_total} todo={len(todo)} n_lhs={n_lhs}", file=sys.stderr)
for i in todo:
    print(i)
PY

n_todo=$(grep -c '^[0-9]' "$to_do_file" || true)
ts=$(date +%Y%m%d_%H%M%S)
jids_file="launch_union_${ts}.jids"
echo "Submitting $n_todo single-sample drivers, ${wave_size}/wave, ${wave_delay}s between waves"
echo "JIDs will be saved to $jids_file"

count=0
n_lhs=$(grep -oE 'n_lhs=[0-9]+' "$to_do_file" | head -1 | cut -d= -f2)
n_lhs=${n_lhs:-5000}

while read -r sidx; do
    [ -z "$sidx" ] && continue
    out=$(sbatch wrap_synthetic_factorial.sh simulate "$mode" "$sidx" "$n_lhs")
    jid=$(echo "$out" | awk '{print $4}')
    echo "$jid $sidx" >> "$jids_file"
    count=$((count + 1))
    if (( count % wave_size == 0 )); then
        echo "  submitted $count/$n_todo (last jid $jid), sleeping ${wave_delay}s"
        sleep "$wave_delay"
    fi
done < "$to_do_file"

rm -f "$to_do_file"

echo
echo "Done submitting $count drivers. JIDs saved to $jids_file"
echo "Monitor:  squeue -u \$USER -o '%.10i %.20j %.8T %.10M' | grep synfac"
echo "When all complete, run:  sbatch wrap_synthetic_factorial.sh plot union"
