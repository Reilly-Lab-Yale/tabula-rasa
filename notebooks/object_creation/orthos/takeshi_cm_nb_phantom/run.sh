#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --job-name=tak_cm_nb
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

# Run python in background so we can trap signals and forward them,
# giving performance_report.__exit__ a chance to write the HTML.
python3 fit.py &
PID=$!

forward_signal() {
    echo "Caught signal $1, forwarding to PID $PID" >&2
    kill -s "$1" "$PID" 2>/dev/null
    wait "$PID"
    EXIT_CODE=$?
    echo "Python exited with code $EXIT_CODE after signal $1" >&2
    exit $EXIT_CODE
}

trap 'forward_signal TERM' SIGTERM
trap 'forward_signal INT'  SIGINT
trap 'forward_signal HUP'  SIGHUP

# Wait for python -- if it exits on its own, capture the code
wait "$PID"
EXIT_CODE=$?
echo "Python exited with code $EXIT_CODE" >&2
exit $EXIT_CODE
