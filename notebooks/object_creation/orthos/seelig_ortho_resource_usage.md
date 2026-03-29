# Seelig Ortho Model Fitting — Resource Usage Summary

**Dataset**: Seelig scMPRA (seelig_scmpra_umiwise.tsv.gz)
**Output**: `seelig_ortho_20260320` (merged from two phases)
**Cluster**: Bouchet (Yale YCRC)
**Account**: prio_skr2
**Generated**: 2026-03-26

---

## Overview

Fitting was split into two independent phases to enable GPU acceleration for
`by_cell_type` models without GPU idle preemption:

| Phase | Model set | Compute | Wall time | Date |
|---|---|---|---|---|
| Phase 1 | `by_cre` | CPU only | 1h 38m | 2026-03-25 |
| Phase 2 | `by_cell_type` | CPU setup + GPU fits | 1h 58m | 2026-03-26 |

Results were merged post-hoc into `seelig_ortho_20260320/`.

---

## Phase 1: by_cre Models (CPU)

### Driver job — 6469126

| Field | Value |
|---|---|
| Job name | seelig_cre |
| State | COMPLETED |
| Partition | priority |
| Node | a1132u07n04 |
| CPUs allocated | 1 |
| Memory requested | 16 GB |
| Submit | 2026-03-25 16:21:52 |
| Start | 2026-03-25 16:21:53 |
| End | 2026-03-25 18:00:13 |
| Wall time | 1:38:20 |
| Time limit | 12:00:00 |
| CPU time used | 18:01 (user: 17:24) |
| CPU efficiency | 18.2% |
| Peak RSS | 2.3 GB |
| Memory efficiency | 8% of 16 GB |
| Disk read | 166.8 MB |
| Disk write | 5.8 MB |
| Exit code | 0:0 |

**Note on CPU efficiency**: Driver job is the Dask scheduler/coordinator; actual
compute ran on the worker node (6469134). Low CPU% on driver is expected.

### CPU worker — 6469134 (1 worker, SLURMCluster)

| Field | Value |
|---|---|
| Job name | seelig_cre_worker |
| State | CANCELLED (by driver at completion) |
| Partition | priority |
| Node | a1130u07n04 |
| CPUs allocated | 1 |
| Memory requested | 60 GB |
| Start | 2026-03-25 16:22:11 |
| End | 2026-03-25 18:00:08 |
| Wall time | 1:37:57 |
| CPU time used | 1:02:59 (user: 56:04) |
| CPU efficiency | 64.3% |
| Peak RSS | 15.0 GB |
| Memory efficiency | 25% of 60 GB |
| Disk read | 38.3 GB |
| Disk write | 189.4 MB |
| Exit code | 0:15 (SIGTERM on cancel — expected) |

**Note**: Disk read of 38.3 GB reflects repeated loading of the data matrix
across TensorZINB fit calls (one per CRE, with `consider_missing=True`).
High disk I/O is the primary bottleneck for this phase; see project notes
on potential Dask-level caching optimisation.

### Phase 1 totals

| Metric | Value |
|---|---|
| Total wall time | ~1h 38m |
| Total CPU-hours consumed | ~1.7 CPU-hours |
| Peak memory (worker) | 15.0 GB |
| Peak disk read | ~38.3 GB |

---

## Phase 2: by_cell_type Models (CPU setup + GPU fits)

Architecture: 1 driver + 4 CPU workers (setup/data prep) + 2 GPU workers
(TensorZINB fits via `device_type="GPU"`).

### Driver job — 6512633

| Field | Value |
|---|---|
| Job name | seelig_ct |
| State | COMPLETED |
| Partition | priority |
| Node | a1132u07n02 |
| CPUs allocated | 1 |
| Memory requested | 16 GB |
| Submit | 2026-03-26 17:21:26 |
| Start | 2026-03-26 17:22:36 |
| End | 2026-03-26 19:20:20 |
| Wall time | 1:57:44 |
| Time limit | 2-00:00:00 |
| CPU time used | 0:49 |
| CPU efficiency | 0.7% |
| Peak RSS | 634 MB |
| Memory efficiency | 3% of 16 GB |
| Disk read | 170.9 MB |
| Disk write | 51.8 MB |
| Exit code | 0:0 |

**Note**: Driver is pure coordinator (Dask scheduler); near-zero CPU is expected.

### CPU workers — 6512677, 6512678, 6512679, 6512680 (4 workers, SLURMCluster)

All workers ran 2026-03-26 17:25–19:20 (~1h 55m). Cancelled by driver at job
completion (exit code 0:15 = SIGTERM, expected).

| Job | Node | CPU efficiency | Peak RSS | Disk read |
|---|---|---|---|---|
| 6512677 | a1130u24n02 | ~29% (33:39 / 1:55:08) | 16.6 GB | 501.3 MB |
| 6512678 | a1130u31n02 | ~0.7% (1:16 / 1:55:07) | 217.8 MB | 136.9 MB |
| 6512679 | a1130u31n02 | ~0.7% (1:16 / 1:55:07) | 365.3 MB | 137.8 MB |
| 6512680 | a1130u31n04 | ~0.8% (1:36 / 1:55:07) | 17.3 GB | 304.1 MB |

**Note**: Worker 6512677 and 6512680 did the bulk of data preparation (high RSS
and disk I/O); workers 6512678/6512679 were largely idle after setup. The
uneven load is expected given Dask's task distribution across the setup phase.

All workers: CPUs allocated = 1, memory requested = 60 GB each.

### GPU workers — 6512773, 6512774 (2 workers, H200, manual sbatch)

Both workers ran on `priority_gpu`, NVIDIA H200 (140.4 GB VRAM).
Cancelled by driver's `finally` block at job completion (expected).

| Field | 6512773 | 6512774 |
|---|---|---|
| Node | a1124u28n01 | a1126u19n01 |
| Start | 2026-03-26 17:35:50 | 2026-03-26 17:40:03 |
| End | 2026-03-26 19:20:18 | 2026-03-26 19:20:18 |
| Wall time | 1:44:28 | 1:40:15 |
| Time limit | 1-00:00:00 | 1-00:00:00 |
| CPU efficiency | 86.5% (1:30:21 / 1:44:28) | 35.0% (35:03 / 1:40:15) |
| Peak CPU RSS | 25.5 GB | 9.1 GB |
| Memory requested | 64 GB | 64 GB |
| **GPU utilization** | **28.1%** | **19.3%** |
| **GPU memory used** | **138.7 GB / 140.4 GB (98.8%)** | **138.7 GB / 140.4 GB (98.8%)** |
| GPUs allocated | 1× H200 | 1× H200 |
| Disk read | 277.1 MB | 274.1 MB |

**Note on GPU utilization**: 19–28% average GPU utilization is consistent with
the ZINB model structure — early-stopping optimizer steps are interspersed with
CPU-side convergence checks, and TF's eager graph execution introduces
synchronisation overhead. Near-100% VRAM utilisation (138.7/140.4 GB) confirms
TensorFlow claimed the full H200 memory for model weights and activations.

Worker 6512773 handled more cell types (higher CPU and GPU utilization);
6512774 had a shorter effective compute window due to later start (waiting for
H200 allocation).

### Phase 2 totals

| Metric | Value |
|---|---|
| Total wall time | ~1h 58m |
| CPU worker CPU-hours (4 workers × ~1:55) | ~7.7 CPU-hours |
| GPU worker CPU-hours (2 workers) | ~2.1 CPU-hours |
| GPU-hours consumed | ~3.4 GPU-hours (H200) |
| Peak CPU memory (worker) | ~25.5 GB (GPU worker 6512773) |
| Peak GPU VRAM | 138.7 GB / 140.4 GB per GPU |
| Peak disk read (CPU worker) | 501 MB |

---

---

## Attempt 15: by_cell_type Refit (2026-03-27) — NB init via tensorzinb v0.0.6

**Supersedes attempt 12 by_cell_type phase.** Refit required after MoM init was
found to give zi_init=0.63 (true ZI ~0.007) for Seelig, causing 5000-step runs
without early-stopping. Fix: MoM fallback in scMPRAforge + `_nb_init()` GPU
support in tensorzinb-plusplus v0.0.6.

The `by_cre` phase (job 6469126, attempt 12) was not rerun — it used NB init
throughout and converged correctly.

### Driver job — 6579440

| Field | Value |
|---|---|
| Job name | seelig_ct |
| State | COMPLETED |
| Partition | priority |
| Node | a1132u07n02 |
| CPUs allocated | 1 |
| Memory requested | 16 GB |
| Submit | 2026-03-27 20:50:47 |
| Start | 2026-03-27 20:51:30 |
| End | 2026-03-27 22:44:28 |
| Wall time | 1:52:58 |
| Time limit | 2-00:00:00 |
| CPU time used | 0:45 |
| CPU efficiency | 0.7% |
| Peak RSS | 1.6 GB |
| Disk read | 170 MB |
| Disk write | 2.8 MB |
| Exit code | 0:0 |

### CPU workers — 6579501–6579504 (4 workers, SLURMCluster)

All workers ran 2026-03-27 20:54–22:44 (~1h 50m). Cancelled by driver at
completion (exit code 0:15 = SIGTERM, expected). Memory requested: 60 GB each.

| Job | Node | CPU efficiency | Peak RSS | Disk read |
|---|---|---|---|---|
| 6579501 | a1130u07n03 | 60.1% (1:06:40 / 1:50:21) | 20.8 GB | 216 MB |
| 6579502 | a1130u09n04 | 1.4% (1:32 / 1:50:21) | 1.0 GB | 505 MB |
| 6579503 | a1130u09n04 | 1.2% (1:17 / 1:50:21) | 358 MB | 136 MB |
| 6579504 | a1130u18n02 | 1.2% (1:21 / 1:50:21) | 13.4 GB | 217 MB |

Worker 6579501 did most of the CPU work (consider_missing setup, design matrix
construction). Workers 6579502/6579503 were largely idle after setup. Worker
6579504 held ~13 GB of data in memory but did minimal CPU compute.

Compared with attempt 12 CPU workers: very similar pattern (one high-RSS worker,
others idle), with similar disk reads.

### GPU workers — 6579584, 6579585 (2 workers, H200, just-in-time sbatch)

Both workers submitted at 20:55:57 via `pre_fit_hook`, ran on `priority_gpu`.
NVIDIA H200 SXM 141 GB VRAM. Cancelled by driver's `finally` block.

| Field | 6579584 | 6579585 |
|---|---|---|
| Node | a1122u11n01 (GPU 6) | a1124u11n01 (GPU 4) |
| Start | 2026-03-27 20:55:57 | 2026-03-27 20:55:57 |
| End | 2026-03-27 22:44:27 | 2026-03-27 22:44:27 |
| Wall time | 1:48:29 | 1:48:29 |
| Time limit | 1-00:00:00 | 1-00:00:00 |
| CPU efficiency | 58.4% (1:03:32 / 1:48:29) | 37.7% (40:59 / 1:48:29) |
| Peak CPU RSS | 23.7 GB | 14.3 GB |
| Memory requested | 64 GB | 64 GB |
| **GPU utilization** | **14.2%** | **16.3%** |
| **GPU memory used** | **138.7 GB / 140.4 GB (98.8%)** | **138.7 GB / 140.4 GB (98.8%)** |
| Disk read | 284 MB | 268 MB |

**GPU utilization note**: 14–16% average (vs 19–28% in attempt 12). Lower
average is consistent with NB init: GPU runs a fast NB-only fit first, then
ZINB optimizer converges quickly from a good starting point rather than grinding
5000 steps from zi=0.63. Near-100% VRAM utilization unchanged — TF claims the
full H200 at import regardless.

### Attempt 15 totals (by_cell_type only)

| Metric | Value |
|---|---|
| Total wall time | 1:52:58 |
| CPU worker CPU-hours (4 × ~1:50) | ~7.3 CPU-hours |
| GPU worker CPU-hours (2 × ~1:48) | ~3.6 CPU-hours |
| GPU-hours consumed | ~3.6 GPU-hours (H200) |
| Peak CPU memory (worker) | 23.7 GB (GPU worker 6579584) |
| Peak GPU VRAM | 138.7 GB / 140.4 GB per GPU |
| Raw stats report | `run_stats_20260327_203757.txt` |

---

## Combined Resource Summary

| Metric | Phase 1 (by_cre) | Phase 2 (by_cell_type) | Total |
|---|---|---|---|
| Wall time | 1:38:20 | 1:57:44 | ~3h 36m |
| CPU-hours | ~1.7 | ~9.8 | ~11.5 |
| GPU-hours (H200) | 0 | ~3.4 | ~3.4 |
| Peak worker RSS | 15.0 GB | 25.5 GB | — |
| Peak VRAM | — | 138.7 GB / 140.4 GB | — |
| Peak disk read | 38.3 GB | 501 MB | ~38.8 GB |

**Hardware**: Bouchet cluster (Yale YCRC). CPU nodes: AMD/Intel compute nodes
(priority partition). GPU nodes: NVIDIA H200 SXM 141 GB (priority_gpu
partition). Interconnect: InfiniBand.

**Software**: TensorZINB with TensorFlow backend, `device_type="GPU"` for
by_cell_type fits. Dask distributed scheduler with heterogeneous worker pool
(CPU-only workers for data prep, GPU workers for TF training).
