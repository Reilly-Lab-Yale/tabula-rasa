#!/usr/bin/env python
"""Median-replicate ROC and PRC curves, one pair per dataset.

Scriptified from the paper-figure cell of all_prc_summary.ipynb so the curves
can be regenerated without stepping through a notebook. Two changes:

  - Wald is dropped from Yin et al. Wald is the only test that reads the
    fitted model rather than the counts, and those orthos were fit under plain
    consider-missing rather than with the MOI correction that dataset's
    canonical fit uses (design fit_mode='cm_phantom'), so its curve there is
    not comparable to the other regimes. Pseudobulk is absent from Yin et al.
    on its own: it needs at least two replicates per group and that design has
    one.
  - The curve shown is chosen the same way as before -- the ground-truth draw
    whose mean auROC is closest to the median across draws -- but the choice
    is asserted rather than assumed.

Needs the simulation objects, hence a Dask client; run it on the cluster.

    python analyses/simulation/activity_prc/plot_curves.py
"""
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge as scm  # noqa: E402

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output" / "paper_figs"
SIM_ROOT = pathlib.Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/simulated")

N_GT_DRAWS = 5
HYPOTHESIS_SET = "hs_all_ct"

DATASETS = {
    "shendure": ("Lalanne et al.", "shendure_5x5_activity"),
    "cohen":    ("Zhao et al.",    "cohen_5x5_activity"),
    "seelig":   ("Yin et al.",     "seelig_5x5_activity"),
}
# Tests to draw per dataset. None means "whatever ran".
TESTS_FOR = {
    "shendure": None,
    "cohen": None,
    "seelig": ["mwu", "ttest", "ks"],
}


def median_gt_draw(slug):
    """The draw whose mean auROC is nearest the median across draws."""
    f = BASE / slug / "output" / f"{slug}_5x5_activity_summary.tsv"
    assert f.is_file(), f"missing summary: {f}"
    d = pd.read_csv(f, sep="\t")
    gt_means = d.groupby("gt_draw")["auroc"].mean()
    assert len(gt_means) == N_GT_DRAWS, (
        f"{slug}: expected {N_GT_DRAWS} ground-truth draws, found {len(gt_means)}")
    pick = int((gt_means - gt_means.median()).abs().idxmin())
    print(f"  {slug}: GT draw {pick} "
          f"(mean auROC {gt_means[pick]:.3f}; draws span "
          f"{gt_means.min():.3f}-{gt_means.max():.3f})")
    return pick


def main():
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=1, threads_per_worker=1, memory_limit="16GB")
    client = Client(cluster)
    OUT.mkdir(parents=True, exist_ok=True)

    for slug, (label, sim_stem) in DATASETS.items():
        gt = median_gt_draw(slug)
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=f"{sim_stem}_gt{gt}", client=client)
        tests = TESTS_FOR[slug]
        for kind in ("ROC", "PRC"):
            sim.median_performance_curve(
                HYPOTHESIS_SET, kind, test_types=tests, include_alpha=True)
            fig = plt.gcf()
            fig.set_size_inches(4.4, 3.6)
            fig.suptitle(f"{label} -- {kind}", fontsize=10)
            for ax in fig.axes:
                if ax.get_legend() is not None:
                    ax.legend(loc="lower right", fontsize=7, framealpha=0.9)
            out = OUT / f"{slug}_median_{kind.lower()}.svg"
            fig.savefig(out, format="svg", bbox_inches="tight")
            plt.close(fig)
            print(f"    wrote {out.name}")

    client.close()
    print("done")


if __name__ == "__main__":
    main()
