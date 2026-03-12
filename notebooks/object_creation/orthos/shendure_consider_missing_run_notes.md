# Shendure Consider-Missing Run Notes

## 2026-03-12 Attempt 1

Driver job:
- Wrapper: [wrap_shend_consider_missing.sh](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/wrap_shend_consider_missing.sh)
- Slurm resources: `2` CPU cores, `24G` RAM, `1-12:10:00` walltime
- Main job id: `1804829`

Dask worker configuration used:
- Defined in [shendure_consider_missing.py](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/shendure_consider_missing.py)
- `SLURMCluster(cores=4, memory="128G", processes=4)`
- `cluster.scale(jobs=4)`
- Effective layout: `16` worker processes total, about `29.8 GiB` memory per worker process

Observed result:
- `ortho_filter` completed
- `consider_missing` modeling ran for roughly `17` hours
- Repeated Dask nanny restarts due to workers exceeding `95%` memory budget
- Final notebook error was `P2PConsistencyError: No active shuffle ... found`
- Root cause was worker loss during shuffle after repeated memory-triggered restarts, not an application-level model exception

Evidence:
- Driver stderr: [slurm-1804829.err](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/slurm-1804829.err)
- Worker logs:
  - [worker_1804850.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804850.out)
  - [worker_1804851.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804851.out)
  - [worker_1804852.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804852.out)
  - [worker_1804853.out](/home/mcn26/project/tabula_rasa/notebooks/object_creation/orthos/worker_1804853.out)

Next configuration to try:
- Reduce per-node Dask process count from `4` to `2`
- Increase cluster size from `4` jobs to `8`
- Keep worker job memory at `128G`
- Add a Dask HTML performance report for memory and shuffle diagnostics
