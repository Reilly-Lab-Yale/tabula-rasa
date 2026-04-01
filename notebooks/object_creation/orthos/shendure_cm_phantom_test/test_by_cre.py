"""
Test: Shendure by_cre CM fit using phantom-zero compression path.

Fits by_cre models through the standard_fit CM path (which now uses
_build_cm_fit_inputs for phantom-zero compression) and compares to
the saved CM ZINB ortho on disk.
"""
import pickle, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-cohen-regen')

from pathlib import Path
from dask.distributed import Client, LocalCluster
import scMPRAforge.core as scm

# ── Paths ──────────────────────────────────────────────────────────────────
DATA_DIR = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data/shendure")
SAVED_ORTHO = DATA_DIR / "shendure_ortho_consider_missing_20260320"

# ── Setup Dask ─────────────────────────────────────────────────────────────
print("Setting up Dask local cluster...", flush=True)
cluster = LocalCluster(n_workers=1, threads_per_worker=2, memory_limit="60GB", processes=False)
client = Client(cluster)
print(f"Dashboard: {client.dashboard_link}", flush=True)

# ── Load data ──────────────────────────────────────────────────────────────
print("Loading Shendure data...", flush=True)
shendure = scm.scMPRA_data.from_tsv(str(DATA_DIR / "shendure_processed.tsv"))
shendure.set_negative_controls(["minP", "noP"])
shendure.set_reference_cell("Pluripotent")
shendure.ortho_filter()
shendure.set_consider_missing(True)
print("Data loaded and filtered.", flush=True)

# ── Load saved ortho for comparison ────────────────────────────────────────
print(f"\nLoading saved CM ortho from {SAVED_ORTHO}...", flush=True)
with open(SAVED_ORTHO / "by_cre.pkl", "rb") as f:
    saved_by_cre = pickle.load(f)

saved_cres = sorted(saved_by_cre.model.keys())
print(f"Saved ortho has {len(saved_cres)} CREs", flush=True)

# Ensure saved weights are numpy
for cre in saved_cres:
    for wk in list(saved_by_cre.model[cre]['weights'].keys()):
        saved_by_cre.model[cre]['weights'][wk] = np.array(
            saved_by_cre.model[cre]['weights'][wk])

# ── Fit by_cre with phantom-zero compressed CM ────────────────────────────
print("\nFitting by_cre with phantom-zero CM path...", flush=True)
ortho_obj = scm.ortho()
t0 = time.time()
ortho_obj.fit_by_cre_models(client=client, dat=shendure)
fit_time = time.time() - t0
print(f"\nFit complete in {fit_time:.1f}s ({fit_time/60:.1f}min)", flush=True)

# ── Gather results ─────────────────────────────────────────────────────────
print("Gathering results...", flush=True)
by_cre = ortho_obj.by_cre
for cre in list(by_cre.model.keys()):
    fut = by_cre.model[cre]
    by_cre.model[cre] = fut.result() if hasattr(fut, 'result') else fut

phantom_cres = sorted(by_cre.model.keys())
print(f"Phantom fit has {len(phantom_cres)} CREs", flush=True)

# ── Compare ────────────────────────────────────────────────────────────────
print("\n" + "="*80, flush=True)
print("COMPARISON: phantom-zero CM vs saved CM ortho (by_cre)", flush=True)
print("="*80, flush=True)

common_cres = sorted(set(saved_cres) & set(phantom_cres))
print(f"Comparing {len(common_cres)} common CREs "
      f"(saved={len(saved_cres)}, phantom={len(phantom_cres)})", flush=True)

llf_diffs = []
mu_diffs = []
pi_diffs = []
theta_diffs = []

for cre in common_cres:
    saved = saved_by_cre.model[cre]
    phantom = by_cre.model[cre]

    # Ensure phantom weights are numpy
    for wk in list(phantom['weights'].keys()):
        phantom['weights'][wk] = np.array(phantom['weights'][wk])

    # LLF comparison
    llf_pct = 100 * (phantom['llf_total'] - saved['llf_total']) / abs(saved['llf_total'])
    llf_diffs.append(llf_pct)

    # x_mu comparison
    s_mu = np.array(saved['weights']['x_mu']).flatten()
    p_mu = np.array(phantom['weights']['x_mu']).flatten()
    denom = np.where(np.abs(s_mu) > 1e-6, s_mu, 1.0)
    mu_pct = np.mean(np.abs((p_mu - s_mu) / denom) * 100)
    mu_diffs.append(mu_pct)

    # x_pi comparison
    if 'x_pi' in saved['weights'] and 'x_pi' in phantom['weights']:
        s_pi = np.array(saved['weights']['x_pi']).flatten()
        p_pi = np.array(phantom['weights']['x_pi']).flatten()
        denom = np.where(np.abs(s_pi) > 1e-6, s_pi, 1.0)
        pi_pct = np.mean(np.abs((p_pi - s_pi) / denom) * 100)
        pi_diffs.append(pi_pct)

    # theta comparison
    s_th = float(np.array(saved['weights']['theta']).flatten()[0])
    p_th = float(np.array(phantom['weights']['theta']).flatten()[0])
    th_pct = 100 * abs(p_th - s_th) / abs(s_th) if abs(s_th) > 1e-10 else 0
    theta_diffs.append(th_pct)

llf_diffs = np.array(llf_diffs)
mu_diffs = np.array(mu_diffs)
pi_diffs = np.array(pi_diffs)
theta_diffs = np.array(theta_diffs)

print(f"\nLLF %diff:   mean={llf_diffs.mean():+.4f}%, "
      f"std={llf_diffs.std():.4f}%, max |%diff|={np.abs(llf_diffs).max():.4f}%", flush=True)
print(f"x_mu |%diff|: mean={mu_diffs.mean():.2f}%, "
      f"median={np.median(mu_diffs):.2f}%, max={mu_diffs.max():.2f}%", flush=True)
if len(pi_diffs) > 0:
    print(f"x_pi |%diff|: mean={pi_diffs.mean():.2f}%, "
          f"median={np.median(pi_diffs):.2f}%, max={pi_diffs.max():.2f}%", flush=True)
print(f"theta |%diff|: mean={theta_diffs.mean():.2f}%, "
      f"median={np.median(theta_diffs):.2f}%, max={theta_diffs.max():.2f}%", flush=True)

# ── ZI comparison (bounds-style) ──────────────────────────────────────────
print("\n" + "="*80, flush=True)
print("BOUNDS-STYLE ZI COMPARISON", flush=True)
print("="*80, flush=True)


def extract_zi(model_weights, design_mats):
    """Extract per-rep mean ZI from fitted weights."""
    Z = design_mats['zi_regressors']
    zi_names = design_mats.get('zi_regressor_names',
                               list(Z.columns) if hasattr(Z, 'columns') else [])
    import scipy.sparse as sp
    Z_arr = Z.toarray() if sp.issparse(Z) else np.array(Z)
    x_pi = np.array(model_weights['x_pi']).flatten()
    linear_zi = Z_arr @ x_pi
    zi_pred = 1 / (1 + np.exp(-linear_zi))
    rep_idx = np.argmax(Z_arr, axis=1)
    rep_labels = []
    for j in rep_idx:
        name = zi_names[j] if j < len(zi_names) else str(j)
        # Extract rep name from "C(rep_id)[xxx]" format
        if '[' in name:
            import re
            m = re.search(r'\[(.*?)\]', name)
            name = m.group(1) if m else name
        rep_labels.append(name)
    df = pd.DataFrame({'rep_id': rep_labels, 'zi': zi_pred})
    return df.groupby('rep_id')['zi'].mean()


# Need design matrices for ZI extraction — use the saved ones
print("Loading saved design matrices for ZI extraction...", flush=True)
saved_design_dir = SAVED_ORTHO / "by_cre_design"
import json
with open(saved_design_dir / "_keys.json") as f:
    design_keys = json.load(f)
all_designs = {}
for i, k in enumerate(design_keys):
    with open(saved_design_dir / f"{i}.pkl", "rb") as f:
        all_designs[k] = pickle.load(f)

phantom_zi_per_cre = {}
saved_zi_per_cre = {}

for cre in common_cres:
    if cre not in all_designs:
        continue
    design = all_designs[cre]
    if 'x_pi' in by_cre.model[cre]['weights']:
        phantom_zi_per_cre[cre] = extract_zi(by_cre.model[cre]['weights'], design)
    if 'x_pi' in saved_by_cre.model[cre]['weights']:
        saved_zi_per_cre[cre] = extract_zi(saved_by_cre.model[cre]['weights'], design)

if phantom_zi_per_cre and saved_zi_per_cre:
    phantom_zi_all = pd.concat(
        [phantom_zi_per_cre[c].rename(c) for c in phantom_zi_per_cre], axis=1)
    saved_zi_all = pd.concat(
        [saved_zi_per_cre[c].rename(c) for c in saved_zi_per_cre], axis=1)

    phantom_by_cre_zi = phantom_zi_all.mean(axis=1)
    saved_by_cre_zi = saved_zi_all.mean(axis=1)

    print("\nPer-rep by_cre_zi comparison:", flush=True)
    print(f"{'rep':<6} {'Saved CM':>12} {'Phantom CM':>12} {'%diff':>12}", flush=True)
    print("-"*46, flush=True)
    for rep in sorted(phantom_by_cre_zi.index):
        sv = saved_by_cre_zi.get(rep, float('nan'))
        ph = phantom_by_cre_zi.get(rep, float('nan'))
        pct = 100 * (ph - sv) / abs(sv) if abs(sv) > 1e-10 else 0
        print(f"{rep:<6} {sv:>12.6f} {ph:>12.6f} {pct:>+11.4f}%", flush=True)

    print(f"\nMean across reps:", flush=True)
    print(f"  Saved CM:    {saved_by_cre_zi.mean():.6f}", flush=True)
    print(f"  Phantom CM:  {phantom_by_cre_zi.mean():.6f}", flush=True)
    pct_mean = 100 * (phantom_by_cre_zi.mean() - saved_by_cre_zi.mean()) / abs(saved_by_cre_zi.mean())
    print(f"  %diff:       {pct_mean:+.4f}%", flush=True)

    # Per-CRE correlation
    common_zi_cres = sorted(set(phantom_zi_per_cre.keys()) & set(saved_zi_per_cre.keys()))
    phantom_cre_means = phantom_zi_all.mean(axis=0)
    saved_cre_means = saved_zi_all.mean(axis=0)
    common_cols = phantom_cre_means.index.intersection(saved_cre_means.index)
    r = np.corrcoef(phantom_cre_means[common_cols], saved_cre_means[common_cols])[0, 1]
    print(f"\nPer-CRE ZI Pearson r: {r:.4f}", flush=True)
    print(f"Per-CRE ZI std: saved={saved_cre_means[common_cols].std():.4f}, "
          f"phantom={phantom_cre_means[common_cols].std():.4f}", flush=True)

# ── Save comparison ────────────────────────────────────────────────────────
import os
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
comparison = {
    'llf_diffs': llf_diffs,
    'mu_diffs': mu_diffs,
    'pi_diffs': pi_diffs,
    'theta_diffs': theta_diffs,
    'fit_time_s': fit_time,
    'n_cres': len(common_cres),
}
if phantom_zi_per_cre and saved_zi_per_cre:
    comparison['phantom_by_cre_zi'] = phantom_by_cre_zi
    comparison['saved_by_cre_zi'] = saved_by_cre_zi

with open(f'{OUT_DIR}/comparison_results.pkl', 'wb') as f:
    pickle.dump(comparison, f)

print(f"\nResults saved to {OUT_DIR}/comparison_results.pkl", flush=True)
print("Done.", flush=True)

client.close()
cluster.close()
