"""Per-fit compute requirements (wall time and peak memory) for the stratified
scMPRA fits.

Parses the SLURM run_stats_*.txt reports emitted for each canonical/
counterfactual stratified fit under analyses/model_fitting/fits/ and produces:

  output/compute_stats.csv        -- intermediate table, one row per fit
  output/compute_scatter.{svg,png}  -- all fits: peak RSS + wall time vs
                                       (CRE, cell-type) groups fit
  output/compute_canonical.{svg,png} -- friendly main-text version: one bar
                                       per source dataset (canonical fit only)

Two regimes are distinguished in the scatter:
  - reporter (obs / obsingle expansions): structural zeros removed via the
    transfection reporter before fitting -> low memory.
  - phantom-zero (cm expansion): a phantom zero per unobserved (cell,
    barcode) pair -> high memory.

All logged fits ran on 4 CPU cores (no GPU); the accelerator axis is not
covered by these logs and is intentionally omitted.
"""
import re
import csv
from pathlib import Path

import matplotlib.pyplot as plt

# Emit editable <text> so the manuscript sync pipeline's global font-size
# scale (and any text substitution) applies, matching the other figures.
plt.rcParams["svg.fonttype"] = "none"

BASE = Path(__file__).parent
FITS = BASE.parent / "model_fitting" / "fits"
OUT = BASE / "output"

REPORTER_EXPANSIONS = {"obs", "obsingle"}  # phantom-zero is "cm"

# Canonical fit per source dataset (see manuscript Methods, model-selection).
# One NB model per dataset. seelig's canonical preset name has no logged fit
# dir (only the zinb variant / non-moib nb exist); fall back to the nearest NB
# fit and flag it.
CANONICAL = {
    "shendure": "shendure_obs_nb_phantom",
    "cohen": "cohen_obsingle_nb_phantom",
    "seelig": "seelig_cm_moib_nb_phantom",
}
CANONICAL_FALLBACK = {"seelig_cm_moib_nb_phantom": "seelig_cm_nb_phantom"}

# code name -> citation label (matches the manuscript's dataset table)
DISPLAY = {
    "shendure": "Lalanne et al.",
    "cohen": "Zhao et al.",
    "seelig": "Yin et al.",
}
# order for the canonical bars (ascending fit size)
CANONICAL_ORDER = ["cohen", "shendure", "seelig"]


def hms_to_seconds(s):
    """Parse sacct elapsed/CPU strings: 'HH:MM:SS', 'H:MM:SS', or 'MM:SS.sss'."""
    s = s.strip()
    if not s:
        return None
    parts = s.split(":")
    parts = [float(p) for p in parts]
    if len(parts) == 3:
        h, m, sec = parts
    elif len(parts) == 2:
        h, m, sec = 0.0, parts[0], parts[1]
    else:
        return parts[0]
    return h * 3600 + m * 60 + sec


def _rss_to_gb(tok):
    m = re.match(r"([\d.]+)([KMG]?)", tok.strip())
    if not m:
        return None
    factor = {"K": 2**-20, "M": 2**-10, "G": 1.0, "": 2**-30}[m.group(2)]
    return float(m.group(1)) * factor


def parse_fit(path):
    txt = path.read_text()
    rec = {"fit": path.parent.name}
    tfit = re.findall(r"Total fit time:\s*([\d.]+)s", txt)
    if not tfit:
        return None  # no successful fit in this report
    rec["total_fit_s"] = float(tfit[-1])
    bycre = re.findall(r"by_cre done in\s*([\d.]+)s", txt)
    rec["by_cre_s"] = float(bycre[-1]) if bycre else None
    byct = re.findall(r"by_cell_type done in\s*([\d.]+)s", txt)
    rec["by_cell_type_s"] = float(byct[-1]) if byct else None
    drop = re.findall(r"Dropped\s+(\d+)\s+of\s+(\d+)\s+\(cell_type, cre_id\) combos", txt)
    if not drop:
        return None
    dropped, total = int(drop[-1][0]), int(drop[-1][1])
    rec["n_groups"] = total - dropped
    # resource fields from the completed .batch sacct record
    peak_rss = avg_rss = total_cpu_s = cpu_reserved_s = elapsed_s = None
    for line in txt.splitlines():
        if re.match(r"^\d+\.batch\|", line) and "COMPLETED" in line:
            f = line.split("|")
            if len(f) > 23 and f[20].strip():
                peak_rss = _rss_to_gb(f[20])       # MaxRSS
                avg_rss = _rss_to_gb(f[21])        # AveRSS
                total_cpu_s = hms_to_seconds(f[16])  # TotalCPU (actual)
                cpu_reserved_s = hms_to_seconds(f[19])  # CPUTime (cores*elapsed)
                elapsed_s = hms_to_seconds(f[14])   # Elapsed
    rec["peak_rss_gb"] = peak_rss
    rec["avg_rss_gb"] = avg_rss
    rec["total_cpu_s"] = total_cpu_s
    rec["cpu_reserved_s"] = cpu_reserved_s
    rec["elapsed_s"] = elapsed_s
    parts = rec["fit"].split("_")
    rec["dataset"] = parts[0]
    rec["model"] = "zinb" if "zinb" in parts else "nb"
    exp = [p for p in parts[1:] if p not in ("nb", "zinb", "phantom")]
    rec["expansion"] = "_".join(exp)
    rec["regime"] = "reporter" if exp[0] in REPORTER_EXPANSIONS else "phantom-zero"
    return rec


def load():
    rows = []
    for d in sorted(FITS.glob("*")):
        if not d.is_dir():
            continue
        chosen = None
        for f in sorted(d.glob("run_stats_*.txt")):
            if "Total fit time" in f.read_text():
                chosen = f  # prefer a report containing a completed fit
        if chosen is None:
            continue
        rec = parse_fit(chosen)
        if rec:
            rows.append(rec)
    return rows


TABLE_COLS = [
    "fit", "dataset", "expansion", "model", "regime", "n_groups",
    "by_cre_s", "by_cell_type_s", "total_fit_s", "elapsed_s",
    "total_cpu_s", "cpu_reserved_s", "avg_rss_gb", "peak_rss_gb",
]


def write_table(rows):
    OUT.mkdir(exist_ok=True)
    path = OUT / "compute_stats.csv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TABLE_COLS)
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["dataset"], x["expansion"], x["model"])):
            w.writerow({k: r.get(k) for k in TABLE_COLS})
    print(f"wrote {path} ({len(rows)} fits)")


# -- detailed scatter over all fits -------------------------------------------
COLORS = {"reporter": "#2c7fb8", "phantom-zero": "#d95f0e"}
MARKERS = {"nb": "o", "zinb": "^"}


def plot_scatter(rows):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6))
    specs = [
        ("peak_rss_gb", "Peak memory (GB)", axes[0]),
        ("total_fit_s", "Wall-clock fit time (min)", axes[1]),
    ]
    for key, ylab, ax in specs:
        for r in rows:
            y = r[key] / 60.0 if key == "total_fit_s" else r[key]
            ax.scatter(
                r["n_groups"], y,
                c=COLORS[r["regime"]], marker=MARKERS[r["model"]],
                s=55, edgecolors="white", linewidths=0.6, zorder=3,
            )
        ax.set_xlabel("(CRE, cell-type) groups fit")
        ax.set_ylabel(ylab)
        if key == "peak_rss_gb":
            ax.set_yscale("log")  # 3-53 GB spans an order of magnitude
        else:
            ax.set_ylim(0, None)  # linear: everything sits under ~90 min
        ax.grid(True, which="major", alpha=0.25, zorder=0)
    axes[0].axhline(32, ls="--", lw=1.0, color="0.4", zorder=1)
    axes[0].text(axes[0].get_xlim()[1], 32, " 32 GB", va="center", ha="left",
                 fontsize=7, color="0.4")

    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="s", ls="", color=COLORS["reporter"],
               label="reporter (obs/obsingle)"),
        Line2D([], [], marker="s", ls="", color=COLORS["phantom-zero"],
               label="phantom-zero (cm)"),
        Line2D([], [], marker="o", ls="", color="0.3", label="NB"),
        Line2D([], [], marker="^", ls="", color="0.3", label="ZINB"),
    ]
    axes[1].legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"compute_scatter.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'compute_scatter.svg'} and .png")


# -- friendly canonical-only bars ---------------------------------------------
BAR_COLOR = "#0072b2"  # Okabe-Ito blue (house palette; global remap is a no-op)


def select_canonical(rows):
    by_fit = {r["fit"]: r for r in rows}
    chosen = []
    for ds in CANONICAL_ORDER:
        want = CANONICAL[ds]
        r = by_fit.get(want)
        if r is None:
            fb = CANONICAL_FALLBACK.get(want)
            r = by_fit.get(fb)
            if r is None:
                print(f"WARNING: no canonical fit for {ds} ({want} / {fb} missing)")
                continue
            print(f"NOTE: canonical fit '{want}' has no log; using fallback "
                  f"'{fb}' for {ds}.")
        chosen.append((ds, r))
    return chosen


def plot_canonical(rows):
    chosen = select_canonical(rows)
    labels = [DISPLAY[ds] for ds, _ in chosen]
    mem = [r["peak_rss_gb"] for _, r in chosen]
    tmin = [r["total_fit_s"] / 60.0 for _, r in chosen]
    x = range(len(chosen))

    fig, (ax_m, ax_t) = plt.subplots(1, 2, figsize=(7.0, 3.2))
    for ax, vals, ylab, unit in [
        (ax_m, mem, "Peak memory (GB)", "GB"),
        (ax_t, tmin, "Fit time (min)", "min"),
    ]:
        ax.bar(x, vals, color=BAR_COLOR, width=0.62, zorder=3)
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, max(vals) * 1.18)
        ax.grid(True, axis="y", alpha=0.25, zorder=0)
        for xi, v in zip(x, vals):
            ax.text(xi, v, f"{v:.0f} {unit}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"compute_canonical.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'compute_canonical.svg'} and .png")


if __name__ == "__main__":
    rows = load()
    for r in sorted(rows, key=lambda x: x["n_groups"]):
        print(f"{r['fit']:32s} groups={r['n_groups']:5d} "
              f"rss={r['peak_rss_gb']:5.1f}GB t={r['total_fit_s']/60:5.1f}min "
              f"[{r['regime']}/{r['model']}]")
    write_table(rows)
    plot_scatter(rows)
    plot_canonical(rows)
