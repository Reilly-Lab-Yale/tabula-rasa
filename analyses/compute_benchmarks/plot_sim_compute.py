"""Simulation compute cost across campaigns.

Aggregates the distributed (main + worker) SLURM run_stats reports under
analyses/simulation/ and produces:

  output/sim_stats.csv            -- intermediate table, one row per campaign
  output/sim_detail.{svg,png}     -- CPU-hours (log) + peak memory per campaign

Each simulation campaign (analysis x dataset) runs a Dask driver ("MAIN") plus
many single-core "WORKER" jobs. Each job block carries a jobstats JSON:
  - top-level "total_time"  -> elapsed wall seconds
  - per-node "total_time"   -> CPU-seconds actually used
  - per-node "cpus"         -> allocated cores
  - per-node "used_memory"  -> peak RSS bytes for that job

design_space/synthetic_factorial is read from its hand-made bucket summary
(its ~24k mixed-run logs do not scale to per-job parsing).
"""
import re
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"

BASE = Path(__file__).parent
SIM = BASE.parent / "simulation"
OUT = BASE / "output"

JOB_HEADER = re.compile(r"^\s*Job\s+(\d+)\s+\[(MAIN|WORKER)\]")

# analysis type -> Okabe-Ito color
ANALYSIS_COLOR = {
    "activity_prc": "#0072b2",
    "activity_power": "#e69f00",
    "activity_calibration": "#009e73",
    "variant_power": "#cc79a7",
    "design_space": "#d55e00",
}


def parse_jobs(txt):
    lines = txt.splitlines()
    idxs = [i for i, ln in enumerate(lines) if JOB_HEADER.match(ln)]
    idxs.append(len(lines))
    for a, b in zip(idxs, idxs[1:]):
        role = JOB_HEADER.match(lines[a]).group(2)
        block = "\n".join(lines[a:b])
        jstart = block.find("jobstats — JSON")
        if jstart == -1:
            continue
        brace = block.find("{", jstart)
        if brace == -1:
            continue
        depth, end = 0, None
        for i in range(brace, len(block)):
            if block[i] == "{":
                depth += 1
            elif block[i] == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        if end is None:
            continue
        try:
            js = json.loads(block[brace:end])
        except json.JSONDecodeError:
            continue
        nodes = js.get("nodes", {})
        yield {
            "role": role,
            "elapsed_s": js.get("total_time", 0) or 0,
            "cpu_s": sum(n.get("total_time", 0) or 0 for n in nodes.values()),
            "cpus": sum(n.get("cpus", 0) or 0 for n in nodes.values()),
            "peak_mem_gb": max((n.get("used_memory", 0) or 0
                                for n in nodes.values()), default=0) / 2**30,
        }


def aggregate_campaign(analysis, dataset, files):
    jobs = []
    for f in files:
        jobs.extend(parse_jobs(f.read_text()))
    if not jobs:
        return None
    workers = [j for j in jobs if j["role"] == "WORKER"]
    mains = [j for j in jobs if j["role"] == "MAIN"]
    compute_jobs = workers if workers else mains  # where the memory lives
    mem = [j["peak_mem_gb"] for j in compute_jobs] or [0]
    actual = sum(j["cpu_s"] for j in jobs) / 3600.0
    reserved = sum(j["cpus"] * j["elapsed_s"] for j in jobs) / 3600.0
    return {
        "analysis": analysis, "dataset": dataset,
        "n_mains": len(mains), "n_workers": len(workers),
        "actual_cpu_h": round(actual, 1), "reserved_cpu_h": round(reserved, 1),
        "eff_pct": round(100 * actual / reserved) if reserved else 0,
        "avg_worker_gb": round(sum(mem) / len(mem), 1),
        "peak_worker_gb": round(max(mem), 1),
        "max_elapsed_h": round(max(j["elapsed_s"] for j in jobs) / 3600.0, 2),
        "n_logs": len(files),
    }


def parse_synthetic_factorial():
    d = SIM / "design_space" / "synthetic_factorial"
    files = list(d.glob("run_stats_*.txt"))
    if not files:
        return None
    txt = files[0].read_text()
    m = re.search(r"TOTAL\s+actual=([\d.]+)\s*cpu-h\s+reserved=([\d.]+)\s*cpu-h"
                  r"\s+overall eff=(\d+)%", txt)
    peak = re.search(r"worker\s+\d+\s+[\d.]+\s+[\d.]+\s+\d+%\s+([\d.]+)\s+([\d.]+)", txt)
    if not m:
        return None
    return {
        "analysis": "design_space", "dataset": "union_sweep",
        "n_mains": None, "n_workers": None,
        "actual_cpu_h": float(m.group(1)), "reserved_cpu_h": float(m.group(2)),
        "eff_pct": int(m.group(3)),
        "avg_worker_gb": float(peak.group(2)) if peak else None,
        "peak_worker_gb": float(peak.group(1)) if peak else None,
        "max_elapsed_h": None, "n_logs": len(files),
    }


def load():
    rows = []
    for adir in sorted(SIM.glob("*")):
        if not adir.is_dir() or adir.name == "design_space":
            continue
        for ds in sorted(adir.glob("*")):
            files = sorted(ds.glob("run_stats_*.txt")) if ds.is_dir() else []
            if not files:
                continue
            rec = aggregate_campaign(adir.name, ds.name, files)
            if rec:
                rows.append(rec)
    sf = parse_synthetic_factorial()
    if sf:
        rows.append(sf)
    return rows


# Datasets that are unpublished collaborator data -> excluded from the paper
# figures but kept in the "with_takeshi" internal variant.
UNPUBLISHED = {"takeshi"}


def in_paper(row):
    return row["dataset"] not in UNPUBLISHED


COLS = ["analysis", "dataset", "in_paper", "n_mains", "n_workers", "actual_cpu_h",
        "reserved_cpu_h", "eff_pct", "avg_worker_gb", "peak_worker_gb",
        "max_elapsed_h", "n_logs"]


def write_table(rows):
    OUT.mkdir(exist_ok=True)
    path = OUT / "sim_stats.csv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in sorted(rows, key=lambda x: -(x["actual_cpu_h"] or 0)):
            w.writerow({k: (in_paper(r) if k == "in_paper" else r.get(k))
                        for k in COLS})
    print(f"wrote {path} ({len(rows)} campaigns)")


def plot_detail(rows, suffix="", note=""):
    rows = sorted(rows, key=lambda x: x["actual_cpu_h"] or 0)  # biggest on top
    labels = [f"{r['analysis'].replace('_', ' ')} / {r['dataset']}" for r in rows]
    colors = [ANALYSIS_COLOR.get(r["analysis"], "0.5") for r in rows]
    y = range(len(rows))

    fig, (ax_c, ax_m) = plt.subplots(
        1, 2, figsize=(10, 5.2), gridspec_kw={"width_ratios": [1.3, 1]})

    # CPU-hours lollipop, log x
    cpu = [r["actual_cpu_h"] for r in rows]
    ax_c.hlines(list(y), 0.5, cpu, color=colors, lw=2, zorder=2)
    ax_c.scatter(cpu, list(y), color=colors, s=45, zorder=3)
    for yi, v in zip(y, cpu):
        ax_c.text(v * 1.15, yi, f"{v:.0f}", va="center", ha="left", fontsize=7)
    ax_c.set_xscale("log")
    ax_c.set_xlim(0.5, max(cpu) * 2.5)
    ax_c.set_xlabel("Actual CPU-hours")
    ax_c.set_yticks(list(y))
    ax_c.set_yticklabels(labels, fontsize=8)
    ax_c.grid(True, axis="x", which="major", alpha=0.25, zorder=0)

    # peak compute-job memory
    mem = [r["peak_worker_gb"] for r in rows]
    ax_m.hlines(list(y), 0, mem, color=colors, lw=2, zorder=2)
    ax_m.scatter(mem, list(y), color=colors, s=45, zorder=3)
    for yi, v in zip(y, mem):
        ax_m.text(v + max(mem) * 0.02, yi, f"{v:.0f}", va="center", ha="left",
                  fontsize=7)
    ax_m.set_xlim(0, max(mem) * 1.18)
    ax_m.set_xlabel("Peak per-process memory (GB)")
    ax_m.set_yticks(list(y))
    ax_m.set_yticklabels([])
    ax_m.grid(True, axis="x", alpha=0.25, zorder=0)

    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=a.replace("_", " "))
               for a, c in ANALYSIS_COLOR.items()]
    ax_c.legend(handles=handles, fontsize=7, frameon=False,
                loc="lower right", title="analysis", title_fontsize=7)
    total = sum(r["actual_cpu_h"] or 0 for r in rows)
    title = (f"Simulation compute, as-run (actual CPU-h, scheduler retries "
             f"included) -- total {total:.0f} CPU-h")
    if note:
        title += f"\n{note}"
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"sim_detail{suffix}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/('sim_detail'+suffix+'.svg')} and .png (total {total:.0f} CPU-h)")


if __name__ == "__main__":
    rows = load()
    write_table(rows)
    paper_rows = [r for r in rows if in_paper(r)]
    plot_detail(rows, suffix="_with_takeshi",
                note="includes unpublished takeshi collaborator campaigns")
    plot_detail(paper_rows, suffix="_no_takeshi",
                note="paper figures only (takeshi excluded)")
