#!/bin/bash
#SBATCH --partition=priority
#SBATCH --account=prio_skr2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=co_obsin_nb
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --exclude=a1132u18n02

module load miniconda
eval "$(conda shell.bash hook)"
conda activate tz

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

wait "$PID"
EXIT_CODE=$?
echo "Python exited with code $EXIT_CODE" >&2
exit $EXIT_CODE
