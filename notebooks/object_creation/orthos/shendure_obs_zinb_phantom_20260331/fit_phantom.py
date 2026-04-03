"""
Refit the Shendure obs ZINB by_cre ortho using phantom-zero weighted compression.

Loads the existing shendure_ortho_20260306 design matrices, builds phantom-zero
compressed inputs for all 208 CREs, fits ZINB with weighted TensorZINB, and saves
results in ortho-compatible format for bounds comparison.

Usage: python fit_phantom.py
"""
import pickle, time, os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
from collections import defaultdict

sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-cohen-regen')
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tensorzinb-plusplus')
from tensorzinb.tensorzinb import TensorZINB

ORTHO_DIR = '/nfs/roberts/project/pi_skr2/shared/tabula_data/shendure/shendure_ortho_20260306'
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Load saved ortho and design matrices ────────────────────────────────────
print("Loading saved ortho...", flush=True)
with open(f'{ORTHO_DIR}/by_cre.pkl', 'rb') as f:
    saved_ortho = pickle.load(f)

print("Loading design matrices...", flush=True)
with open(f'{ORTHO_DIR}/by_cre_design.pkl', 'rb') as f:
    all_designs = pickle.load(f)

cre_names = sorted(saved_ortho.model.keys())
print(f"Found {len(cre_names)} CREs", flush=True)

# ── Build compressed inputs and fit each CRE ────────────────────────────────

def compress_cre(design):
    """Build phantom-zero compressed inputs from a design matrix dict."""
    nb_csr = sp.csr_matrix(design['nb_regressors'])
    zi_csr = sp.csr_matrix(design['zi_regressors'])
    y_full = design['regressand']['umis_mpra_bc'].values
    n = len(y_full)
    n_nb_cols = nb_csr.shape[1]
    n_zi_cols = zi_csr.shape[1]

    # Identify group membership for each row
    nb_col_indices = np.zeros(n, dtype=np.int32)
    zi_col_indices = np.zeros(n, dtype=np.int32)
    for i in range(n):
        s, e = nb_csr.indptr[i], nb_csr.indptr[i+1]
        tc = nb_csr.indices[s:e]; tc = tc[tc > 0]
        if len(tc): nb_col_indices[i] = tc[0]
        s, e = zi_csr.indptr[i], zi_csr.indptr[i+1]
        zi_col_indices[i] = zi_csr.indices[s]

    groups = defaultdict(lambda: {'total': 0, 'nonzero_y': []})
    for i in range(n):
        key = (nb_col_indices[i], zi_col_indices[i])
        groups[key]['total'] += 1
        if y_full[i] > 0:
            groups[key]['nonzero_y'].append(y_full[i])

    comp_y, comp_w, comp_nb_rows, comp_zi_rows = [], [], [], []
    row_idx = 0
    for (nb_col, zi_col), g in sorted(groups.items()):
        for yval in g['nonzero_y']:
            comp_y.append(yval); comp_w.append(1.0)
            comp_nb_rows.append((row_idx, nb_col))
            comp_zi_rows.append((row_idx, zi_col))
            row_idx += 1
        nz = g['total'] - len(g['nonzero_y'])
        if nz > 0:
            comp_y.append(0); comp_w.append(float(nz))
            comp_nb_rows.append((row_idx, nb_col))
            comp_zi_rows.append((row_idx, zi_col))
            row_idx += 1

    n_c = row_idx
    comp_y = np.array(comp_y, dtype=np.int64)
    comp_w = np.array(comp_w, dtype=np.float64)

    nb_ri, nb_ci, nb_d = [], [], []
    for (r, c) in comp_nb_rows:
        nb_ri.append(r); nb_ci.append(0); nb_d.append(1.0)
        if c > 0: nb_ri.append(r); nb_ci.append(c); nb_d.append(1.0)
    comp_nb = sp.csr_matrix((nb_d, (nb_ri, nb_ci)), shape=(n_c, n_nb_cols))
    comp_zi = sp.csr_matrix(
        ([1.0]*len(comp_zi_rows),
         ([r for r,c in comp_zi_rows], [c for r,c in comp_zi_rows])),
        shape=(n_c, n_zi_cols)
    )

    assert comp_w.sum() == n, f"Weight sum {comp_w.sum()} != {n}"
    return comp_y, comp_w, comp_nb, comp_zi, n, n_c


def extract_zi(model_weights, design):
    """Extract per-rep mean ZI from fitted weights, mirroring _extract_zi in core.py."""
    Z = design['zi_regressors']
    zi_names = list(Z.columns) if hasattr(Z, 'columns') else design.get('zi_regressor_names', [])
    Z_arr = np.array(Z)
    x_pi = np.array(model_weights['x_pi']).flatten()
    linear_zi = Z_arr @ x_pi
    zi_pred = 1 / (1 + np.exp(-linear_zi))

    # Build rep labels from one-hot
    rep_labels = []
    for i in range(Z_arr.shape[0]):
        idx = np.argmax(Z_arr[i])
        rep_labels.append(zi_names[idx].replace('C(rep_id)[', '').replace(']', ''))

    df = pd.DataFrame({'rep_id': rep_labels, 'zi': zi_pred})
    return df.groupby('rep_id')['zi'].mean()


# ── Main fitting loop ───────────────────────────────────────────────────────
phantom_results = {}  # CRE -> result dict (same format as ortho.model[CRE])
phantom_zi_per_cre = {}  # CRE -> Series(rep_id -> mean_zi)
saved_zi_per_cre = {}

total_t0 = time.time()

for i, cre in enumerate(cre_names):
    design = all_designs[cre]
    saved = saved_ortho.model[cre]

    # Convert saved weights to numpy
    for wk in list(saved['weights'].keys()):
        saved['weights'][wk] = np.array(saved['weights'][wk])

    comp_y, comp_w, comp_nb, comp_zi, n_full, n_comp = compress_cre(design)

    t0 = time.time()
    zinbo = TensorZINB(
        endog=comp_y, exog=comp_nb, exog_infl=comp_zi,
        sample_weight=comp_w, nb_only=False
    )
    res = zinbo.fit(return_history=True, init_method="nb", device_type="CPU")
    fit_time = time.time() - t0

    # Label weights with regressor names (matching ortho format)
    nb_names = list(design['nb_regressors'].columns)
    zi_names = list(design['zi_regressors'].columns)
    res['weights']['x_mu'] = pd.Series(res['weights']['x_mu'].flatten(), index=nb_names)
    res['weights']['x_pi'] = pd.Series(res['weights']['x_pi'].flatten(), index=zi_names)

    phantom_results[cre] = res

    # Extract ZI for bounds comparison
    phantom_zi_per_cre[cre] = extract_zi(res['weights'], design)
    saved_zi_per_cre[cre] = extract_zi(saved['weights'], design)

    # LLF comparison
    llf_pct = 100 * (res['llf_total'] - saved['llf_total']) / abs(saved['llf_total'])

    if (i + 1) % 10 == 0 or i == 0:
        print(f"[{i+1:3d}/{len(cre_names)}] {cre:<25} "
              f"n={n_full:>6,}→{n_comp:>4,} ({n_full/n_comp:.0f}x)  "
              f"LLF%={llf_pct:+.4f}%  "
              f"epochs={res['epochs']}  t={fit_time:.1f}s", flush=True)
    del zinbo

total_time = time.time() - total_t0
print(f"\nAll {len(cre_names)} CREs fit in {total_time:.1f}s ({total_time/60:.1f}min)", flush=True)

# ── Save phantom ortho (by_cre only) ────────────────────────────────────────
# We save as a simple dict since we don't need the full experiment_model class
# for bounds comparison. But we match the format.
print("\nSaving results...", flush=True)

# Save per-CRE results
with open(f'{OUT_DIR}/phantom_by_cre_results.pkl', 'wb') as f:
    pickle.dump(phantom_results, f)

# ── Compute and compare bounds-style ZI ─────────────────────────────────────
print("\n" + "="*80, flush=True)
print("BOUNDS-STYLE ZI COMPARISON", flush=True)
print("="*80, flush=True)

# Average ZI across all CREs per rep (same as Bounds.from_ortho)
phantom_zi_dfs = []
saved_zi_dfs = []
for cre in cre_names:
    pzi = phantom_zi_per_cre[cre].rename(cre)
    szi = saved_zi_per_cre[cre].rename(cre)
    phantom_zi_dfs.append(pzi)
    saved_zi_dfs.append(szi)

phantom_zi_all = pd.concat(phantom_zi_dfs, axis=1)  # rows=rep_id, cols=CRE
saved_zi_all = pd.concat(saved_zi_dfs, axis=1)

phantom_by_cre_zi = phantom_zi_all.mean(axis=1)
saved_by_cre_zi = saved_zi_all.mean(axis=1)

# Load existing SHENDURE_BOUNDS for reference
from scMPRAforge.core import SHENDURE_BOUNDS

print("\nPer-rep by_cre_zi comparison:", flush=True)
print(f"{'rep':<6} {'SHENDURE_BOUNDS':>16} {'Saved ortho':>12} {'Phantom':>12} {'%diff(phantom vs saved)':>24}", flush=True)
print("-"*75, flush=True)

for rep in sorted(phantom_by_cre_zi.index):
    sb = SHENDURE_BOUNDS.by_cre_zi.get(rep, float('nan'))
    sv = saved_by_cre_zi.get(rep, float('nan'))
    ph = phantom_by_cre_zi.get(rep, float('nan'))
    pct = 100 * (ph - sv) / abs(sv) if abs(sv) > 1e-10 else 0
    print(f"{rep:<6} {sb:>16.6f} {sv:>12.6f} {ph:>12.6f} {pct:>+23.4f}%", flush=True)

print(f"\nMean across reps:", flush=True)
print(f"  SHENDURE_BOUNDS: {SHENDURE_BOUNDS.by_cre_zi.mean():.6f}", flush=True)
print(f"  Saved ortho:     {saved_by_cre_zi.mean():.6f}", flush=True)
print(f"  Phantom:         {phantom_by_cre_zi.mean():.6f}", flush=True)
print(f"  %diff:           {100*(phantom_by_cre_zi.mean() - saved_by_cre_zi.mean())/abs(saved_by_cre_zi.mean()):+.4f}%", flush=True)

# Also show per-CRE ZI spread
print(f"\nPer-CRE ZI spread (averaged across reps):", flush=True)
phantom_cre_means = phantom_zi_all.mean(axis=0)  # mean across reps for each CRE
saved_cre_means = saved_zi_all.mean(axis=0)
pct_diffs = 100 * (phantom_cre_means - saved_cre_means) / saved_cre_means.where(saved_cre_means.abs() > 1e-10, 1.0)

print(f"  Mean |%diff| per CRE:   {pct_diffs.abs().mean():.2f}%", flush=True)
print(f"  Median |%diff| per CRE: {pct_diffs.abs().median():.2f}%", flush=True)
print(f"  Max |%diff| per CRE:    {pct_diffs.abs().max():.2f}% ({pct_diffs.abs().idxmax()})", flush=True)
print(f"  Std of %diff:           {pct_diffs.std():.2f}%", flush=True)

# Save comparison
comparison = {
    'phantom_by_cre_zi': phantom_by_cre_zi,
    'saved_by_cre_zi': saved_by_cre_zi,
    'shendure_bounds_by_cre_zi': SHENDURE_BOUNDS.by_cre_zi,
    'phantom_zi_all': phantom_zi_all,
    'saved_zi_all': saved_zi_all,
    'per_cre_pct_diffs': pct_diffs,
    'total_fit_time_s': total_time,
}
with open(f'{OUT_DIR}/bounds_comparison.pkl', 'wb') as f:
    pickle.dump(comparison, f)

print(f"\nResults saved to {OUT_DIR}/", flush=True)
print("Done.", flush=True)
