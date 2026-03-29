# Seelig Ortho Run Graveyard

## 2026-03-24 Attempt 1

Driver job:
- Wrapper: `wrap_fit_seelig_ortho.sh`
- Slurm resources: `1` CPU core, `16G` RAM, `2-00:00:00` walltime, `priority` partition, `prio_skr2` account
- Main job id: `6438142`
- Started: `2026-03-24 18:06`, failed: `18:06` (~2 seconds in)

Observed result:
- Immediate failure with exit code `1:0`

Root cause:
- Wrapper had `set -euo pipefail` and `conda activate tz`
- `conda activate` fails in a non-interactive shell, causing immediate exit under `pipefail`

Fix applied:
- Removed `set -euo pipefail`
- Added `module reset` before `module load miniconda`
- Switched from `python` to `ipython fit_seelig_ortho.py` (matches working shendure pattern)

---

## 2026-03-24 Attempt 2

Driver job:
- Wrapper: `wrap_fit_seelig_ortho.sh`
- Main job id: `6438147`
- Started: `2026-03-24 18:06`, failed: `18:07` (~3 seconds in)

Observed result:
- Shell printed `EXITING SHELL` (exit 0) but no model was created
- ipython error: `[TerminalIPythonApp] WARNING | File 'fit_seelig_ortho.py' doesn't exist`

Root cause:
- Wrapper used `cd "$(dirname "$0")"` to find the script directory
- Under sbatch, `$0` resolves to a temp copy of the script in `/var/spool/slurmd/...`, not the submission directory
- ipython launched from the wrong directory and couldn't find `fit_seelig_ortho.py`

Fix applied:
- Replaced `cd "$(dirname "$0")"` with the absolute path:
  `cd /nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/notebooks/object_creation/orthos`

---

## 2026-03-24–25 Attempt 3

Driver job:
- Wrapper: `wrap_fit_seelig_ortho.sh`
- Slurm resources: `1` CPU core, `16G` RAM, `2-00:00:00` walltime
- Script: `fit_seelig_ortho.py`
- Main job id: `6438149`
- Workers: `6438151`–`6438158`
- Started: `2026-03-24 18:07`, driver cancelled manually: `2026-03-25 11:57` (~17:49 elapsed)

Dask worker configuration:
- `SLURMCluster(cores=1, memory="64G", processes=1)`
- `cluster.scale(jobs=8)`
- Worker walltime: `--time=12:00:00`
- Effective layout: 8 single-threaded worker processes, 64 GiB each

Data / filtering:
- `ortho_filter` removed 1 (cell_type, cre_id) combo involving `'reference'`
- Dropped 10 of 2687 (cell_type, cre_id) combos with fewer than 3 nonzero entries
- `consider_missing=True` (no transfection reporter in seelig)

Observed result:
- By-CRE models (workers 6438151–6438155, 6438158): completed successfully within ~17 minutes of start
- By-cell-type models (workers 6438156, 6438157): MoM initialization completed (~19:16), then entered TF training
- All 8 workers hit the 12-hour wall at `2026-03-25 06:08` and were `TIMEOUT`-cancelled
- By-cell-type TF training was still in progress at 12 hours — did not converge
- Main driver job left stranded with 0 workers; manually cancelled

Root cause:
- Worker walltime of 12h is insufficient for by-cell-type model TF training on CPU
- seelig has 2 cell types (HepG2 + 1 other) → 2 large models, each >12h on CPU
- By comparison, shendure by-cell-type models were estimated >8h; seelig appears similar or worse

Fix applied:
- Added GPU support via Dask resource tokens: `criss_cross(gpu=True)` routes cell-type fits to `resources={"GPU": 1}` workers
- Added 2 manual GPU workers (`priority_gpu`, RTX 5000 Ada) via `sbatch dask-worker <scheduler_addr> --resources GPU=1`

---

## 2026-03-25 Attempt 4

Driver job:
- Main job id: `6461080`; GPU workers: `6461097`, `6461098`
- Started: ~12:13, failed: ~12:13 (<1 min)

Observed result:
- `EOFError: Ran out of input` loading `by_cell_type.pkl`

Root cause:
- Stale partial model directory from attempt 3 left a 0-byte `by_cell_type.pkl`
- Workers were immediately cancelled when driver exited

Fix applied:
- Deleted `seelig_ortho_20260320/` directory

---

## 2026-03-25 Attempt 5

Driver job:
- Main job id: `6461108`; GPU workers: `6461122`, `6461123`
- GPU workers connected ~12:15, cancelled ~13:10 (~55 min idle)

Observed result:
- GPU workers connected and idle for 55 min then preempted by `priority_gpu` cluster policy
- CRE models were still running when GPU workers were cancelled

Root cause:
- `priority_gpu` cancels jobs not using GPU after ~1 hour
- Setup (consider_missing, ~50 min) + partial CRE training kept GPU workers idle past the threshold

---

## 2026-03-25 Attempt 6

Driver job:
- Main job id: `6462613`; GPU workers: `6462623`, `6462624`

Observed result:
- `ModuleNotFoundError: No module named 'scMPRAforge'` on GPU workers when tasks were dispatched

Root cause:
- Dask nanny forks a new subprocess to run the worker; does not inherit shell `cwd`
- `cd {work_dir}` in the sbatch script does not make local packages importable in the forked worker

Fix applied:
- `export PYTHONPATH={work_dir}:$PYTHONPATH` in the GPU sbatch script (inherited by nanny subprocess)

---

## 2026-03-25 Attempt 7

Driver job:
- Main job id: `6463663`; GPU workers: `6463674`, `6463675`
- GPU workers connected ~14:13, preempted ~15:15 (~62 min idle)

Observed result:
- PYTHONPATH fix worked; workers connected and registered successfully
- GPU workers idle for 62 min then preempted by `priority_gpu` cluster policy
- Same root cause as attempt 5: CRE setup+training exhausts idle GPU budget before cell-type fits start

Root cause:
- Total time before GPU tasks are dispatched: ~50 min setup + ~17 min CRE models = ~67 min > ~60 min preemption threshold
- GPU workers submitted at job start, sit idle during entire CRE phase

Fix applied:
- Split into two independent jobs: `fit_seelig_ortho_by_cell_type.py` (GPU) and `fit_seelig_ortho_by_cre.py` (CPU-only)
- GPU job runs `fit_by_cell_type_models` only — no CRE models, so GPU tasks start ~50 min after setup (safely under preemption threshold)
- CPU job runs `fit_by_cre_models` only, in parallel
- Post-hoc merge script combines both saved objects into final `seelig_ortho_20260320`

---

## 2026-03-26 Attempt 8

Driver job:
- Script: `fit_seelig_ortho_by_cell_type.py` (by-cell-type only, no CRE)
- Main job id: `6469040`; GPU workers: `6469045`, `6469046`
- Failed: ~10 min in

Observed result:
- `WorkerStartTimeoutError: Only 4/6 workers arrived after 600`
- GPU workers stuck PENDING (Resources) — never started

Root cause:
- GPU partition resource contention; workers didn't start within 10 min timeout

---

## 2026-03-26 Attempt 9

Driver job:
- Main job id: `6472118`; GPU workers: `6472125`, `6472126`
- Timeout increased to 1200s (20 min)

Observed result:
- `WorkerStartTimeoutError: Only 4/6 workers arrived after 1200`
- Same resource contention; GPU workers still PENDING after 20 min

---

## 2026-03-26 Attempt 10

Driver job:
- Main job id: `6494756`; GPU workers: `6494807`, `6494808` (H200 partition)
- GPU workers connected successfully (~20 min)

Observed result:
- GPU workers connected and idle ~1:22 then CANCELLED at 13:00:22
- Same idle-GPU preemption, now on H200 partition
- Driver remained running but GPU workers gone; fits would never complete

Root cause:
- H200 partition has same idle-GPU preemption policy as RTX 5000 Ada
- Setup (consider_missing, ~50 min) + MoM init left GPU workers idle >60 min before any fit task dispatched

Fix applied:
- Added `pre_fit_hook` parameter to `standard_fit` and `fit_by_cell_type_models` in `core.py`
- Hook fires after `dask.distributed.wait()` on setup futures, before `_tensorzinb_fit` submission
- Script starts CPU-only, then hook sbatches GPU workers and waits for connect just before fits
- GPU workers now receive tasks within seconds of starting, not minutes

---

## 2026-03-26 Attempt 11

Driver job:
- Main job id: `6505567`; GPU workers: `6505819`, `6505820` (H200 partition)
- GPU workers connected successfully, fits submitted
- Cancelled manually after ~1h30m

Observed result:
- GPU workers connected and fits dispatched (pre_fit_hook working)
- 141 GiB GPU memory allocated on both H200s (TF eager allocation)
- nvidia-smi showed 0% GPU utilization sustained over 1 min of 5s polling
- TF was running all compute on CPU despite GPU workers

Root cause:
- `TensorZINB.fit()` defaults to `device_type="CPU"` — GPU is never requested
- TF allocates all available GPU memory at import time regardless of device placement
- `fit_resources={"GPU": 1}` only controls Dask task routing, not TF device placement

Fix applied:
- Added `use_gpu=False` parameter to `_tensorzinb_fit`
- `standard_fit` derives `use_gpu = "GPU" in fit_resources` and passes it through `client.submit`
- `_tensorzinb_fit` sets `device_type = "GPU" if use_gpu else "CPU"` and passes to all `zinbo.fit()` calls

---

## 2026-03-26 Attempt 12 — SUCCESS (superseded)

Driver job: `6512633`; GPU workers: `6512773`, `6512774` (H200); CPU workers: `6512677`–`6512680`

Result:
- Setup complete in ~8 min; GPU workers submitted and connected
- GPU utilization confirmed: 28% (6512773) and 19% (6512774), 138.7/140.4 GB VRAM
- Fits completed in ~1h 58m total wall time
- Output saved to `seelig_ortho_20260320_by_cell_type/`
- Merged with `seelig_ortho_20260320_by_cre/` (job 6469126) into `seelig_ortho_20260320/`
- Resource usage documented in `seelig_ortho_resource_usage.md`

Superseded by attempt 13: MoM initialization found to be poor for low-count datasets.
Analysis showed by_cell_type ZI init = 0.63, true ZI ~0.007–0.018; model hit max steps (5000)
and never early-stopped. by_cre (NB init) converged in ~100 steps to ZI ~0.006–0.007.

---

## 2026-03-27 Attempt 13 — CANCELLED

Driver job: `6578041`; CPU workers: `6578094`–`6578097`; GPU workers: `6578148`, `6578149`

MoM fallback fired correctly (P(X=0|NB)=0.292 for reference, 0.296 for K562 — both > 0.05).
Both cell types correctly fell back to NB init. However GPU utilization remained 0% throughout.

Root cause:
- `TensorZINB._nb_init()` internally creates a new TensorZINB with `nb_only=True` and calls
  `nb_mod.fit(init_method="poi")` without passing `device_type` — always runs on CPU
- Workers were pegged at ~98% CPU doing the NB-only fit before TF optimizer ever started
- In attempt 12 (MoM init), `_nb_init()` was never called — weights passed directly via
  `init_weights`, jumping straight to GPU optimizer. That's why GPU hit 28% immediately.

Fix required:
- Patch `tensorzinb` package: pass `device_type` through `_nb_init()` to `nb_mod.fit()`
- Alternatively: in scMPRAforge MoM fallback, pass precomputed MoM NB weights but override
  ZI to logit(0.01), skipping `_nb_init()` entirely while still getting correct ZI start

Cancelled after ~30 min. Logs swept. `seelig_ortho_20260320_by_cell_type/` not created.

---

## 2026-03-27 Attempt 14 — FAILED (immediate)

Driver job: `6579104`

Observed result:
- `EOFError: Ran out of input` loading `by_cell_type.pkl`

Root cause:
- Attempt 13 cancellation left a 0-byte `by_cell_type.pkl` and 4-byte `by_cre.pkl` in a
  newly-created `seelig_ortho_20260320_by_cell_type/` directory (written at job startup)
- Attempt 14 found the stale dir and tried to load it

Fix: deleted stale `seelig_ortho_20260320_by_cell_type/`. Resubmitted as attempt 15.

---

## 2026-03-27 Attempt 15 — SUCCESS

Driver job: `6579440`; CPU workers: `6579501`–`6579504`; GPU workers: `6579584`, `6579585` (H200)

- Started: ~20:55, completed: ~22:44 (~1h50m total)
- GPU utilization: confirmed active (86%/41% at peak), 141 GB VRAM allocated
- MoM fallback fired correctly: P(X=0|NB)=0.292 (reference), 0.296 (K562) — both fell back to NB init
- GPU phase completed successfully; workers spent ~45 min in CPU serialization before returning results
- `by_cell_type.pkl`: 2.7 MB; `by_cell_type_parameters.pkl`: 424 KB

Post-fit:
- Merged with `seelig_ortho_20260320_by_cre` (job `6469126`, by_cre.pkl=4.0 MB) into `seelig_ortho_20260320/`
- QC job `6583556` submitted from `seelig_ortho_20260320_qc/wrap_seelig_ortho_qc.sh`

Resource usage: `seelig_ortho_resource_usage.md` § "Attempt 15" (raw stats: `run_stats_20260327_203757.txt`)

Key fix stack for this attempt:
- tensorzinb-plusplus v0.0.6: `_nb_init()` now passes `device_name` through → NB init runs on GPU
- scMPRAforge MoM fallback: `_mom_from_training_data()` returns None when median P(X=0|NB) > 0.05
- `_tensorzinb_fit()`: handles `init_vals is None` by falling back to `init_method="nb"`
