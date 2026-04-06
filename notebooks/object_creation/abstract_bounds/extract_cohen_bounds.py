"""
Extract Bounds presets from all 4 cohen phantom orthos.
Saves .tgz presets and transfection/library model plots as SVGs.
"""
import sys
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa')

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dask.distributed import Client, LocalCluster
import scMPRAforge.core as scm
from pathlib import Path

DATA_DIR = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/cohen")
PRESET_DIR = Path("/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/scMPRAforge/presets")
OUT_DIR = Path(__file__).parent / "output" / "cohen"

ORTHOS = [
    "cohen_obs_zinb_phantom_20260401",
    "cohen_obs_nb_phantom_20260401",
    "cohen_cm_zinb_phantom_20260401",
    "cohen_cm_nb_phantom_20260401",
]

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cluster = LocalCluster(n_workers=1, threads_per_worker=1,
                           processes=False, memory_limit="8GB")
    client = Client(cluster)
    print(f"Dashboard: {client.dashboard_link}", flush=True)

    for name in ORTHOS:
        print(f"\n{'='*60}", flush=True)
        print(f"Processing {name}", flush=True)
        print(f"{'='*60}", flush=True)

        ortho = scm.ortho.load(client, str(DATA_DIR), name)
        print(f"  by_cre: {len(ortho.by_cre.model)} models", flush=True)
        print(f"  by_cell_type: {len(ortho.by_cell_type.model)} models", flush=True)

        bounds = scm.Bounds.from_ortho(ortho)
        print(f"  min_mpra_umi: {bounds.min_mpra_umi:.6f}", flush=True)
        print(f"  max_mpra_umi: {bounds.max_mpra_umi:.6f}", flush=True)
        print(f"  by_cell_type_theta: {bounds.by_cell_type_theta:.6f}", flush=True)
        print(f"  by_cre_theta: {bounds.by_cre_theta:.6f}", flush=True)
        print(f"  reference_activity: {bounds.reference_activity:.6f}", flush=True)
        print(f"  rep_ids: {bounds.rep_ids}", flush=True)
        print(f"  consider_missing: {bounds.consider_missing}", flush=True)
        print(f"  nb_only: {bounds.nb_only}", flush=True)
        print(f"  n_negative_controls: {bounds.n_negative_controls}", flush=True)
        print(f"  zi type: {type(bounds.zi).__name__}, is None: {bounds.zi is None}", flush=True)

        # Short name for files (strip date suffix)
        short = name.replace("_20260401", "")

        # Save preset
        tgz_path = PRESET_DIR / f"{short}.tgz"
        bounds.to_tgz(str(tgz_path))
        print(f"  Saved preset: {tgz_path}", flush=True)

        # Transfection model plot
        fig, ax = plt.subplots(figsize=(8, 5))
        plt.sca(ax)
        bounds.plot_transfection()
        ax.set_title(f"{short} -- transfection model")
        fig.tight_layout()
        svg_path = OUT_DIR / f"{short}_transfection.svg"
        fig.savefig(str(svg_path))
        plt.close(fig)
        print(f"  Saved: {svg_path}", flush=True)

        # Library model plot
        fig, ax = plt.subplots(figsize=(8, 5))
        plt.sca(ax)
        bounds.library_model.plot()
        ax.set_title(f"{short} -- library model")
        fig.tight_layout()
        svg_path = OUT_DIR / f"{short}_library.svg"
        fig.savefig(str(svg_path))
        plt.close(fig)
        print(f"  Saved: {svg_path}", flush=True)

    client.close()
    cluster.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
