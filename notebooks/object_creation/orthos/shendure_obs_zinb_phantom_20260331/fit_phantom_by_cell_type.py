"""
Refit the Shendure obs ZINB by_cell_type ortho using phantom-zero weighted compression.

Loads the existing shendure_ortho_20260306 by_cell_type design matrices, builds
phantom-zero compressed inputs for all 10 cell types, fits ZINB with weighted
TensorZINB, and saves results for bounds comparison.

Usage: python fit_phantom_by_cell_type.py
"""
import pickle, time, os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-cohen-regen')
sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tensorzinb-plusplus')
from tensorzinb.tensorzinb import TensorZINB

ORTHO_DIR = '/nfs/roberts/project/pi_skr2/shared/tabula_data/shendure/shendure_ortho_20260306'
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Load saved ortho and design matrices ────────────────────────────────────
print("Loading saved by_cell_type ortho...", flush=True)
with open(f'{ORTHO_DIR}/by_cell_type.pkl', 'rb') as f:
    saved_ortho = pickle.load(f)

print("Loading by_cell_type design matrices...", flush=True)
with open(f'{ORTHO_DIR}/by_cell_type_design.pkl', 'rb') as f:
    all_designs = pickle.load(f)

cell_types = sorted(saved_ortho.model.keys())
print(f"Found {len(cell_types)} cell types", flush=True)


# ── Vectorized compression ──────────────────────────────────────────────────

def compress(design):
    """Build phantom-zero compressed inputs from a design matrix dict.

    Uses vectorized numpy ops for group assignment instead of row-by-row Python.
    """
    nb_csr = sp.csr_matrix(design['nb_regressors'])
    zi_csr = sp.csr_matrix(design['zi_regressors'])
    y_full = design['regressand']['umis_mpra_bc'].values
    n = len(y_full)
    n_nb = nb_csr.shape[1]
    n_zi = zi_csr.shape[1]

    # ── Assign each row to its (nb_col, zi_col) group ───────────────────
    # NB: treatment contrast — each row has intercept (col 0) + at most one
    # treatment column. The "group" nb_col is the treatment column index,
    # or 0 if reference level (intercept only).
    nb_col = np.zeros(n, dtype=np.int32)
    # For each row, find non-intercept nonzero column
    for i in range(n):
        s, e = nb_csr.indptr[i], nb_csr.indptr[i+1]
        cols = nb_csr.indices[s:e]
        treatment = cols[cols > 0]
        if len(treatment):
            nb_col[i] = treatment[0]

    # ZI: one-hot rep indicators — each row has exactly one nonzero entry
    zi_col = np.zeros(n, dtype=np.int32)
    for i in range(n):
        s, e = zi_csr.indptr[i], zi_csr.indptr[i+1]
        zi_col[i] = zi_csr.indices[s]

    # ── Group and compress ──────────────────────────────────────────────
    # Encode group as single int for fast grouping
    group_key = nb_col.astype(np.int64) * n_zi + zi_col.astype(np.int64)
    unique_groups = np.unique(group_key)

    comp_y_list = []
    comp_w_list = []
    comp_nb_col_list = []
    comp_zi_col_list = []

    for gk in unique_groups:
        mask = group_key == gk
        g_nb = int(gk // n_zi)
        g_zi = int(gk % n_zi)
        g_y = y_full[mask]

        # Nonzero observations: individual rows, weight=1
        nz_mask = g_y > 0
        nz_y = g_y[nz_mask]
        if len(nz_y):
            comp_y_list.append(nz_y)
            comp_w_list.append(np.ones(len(nz_y)))
            comp_nb_col_list.extend([g_nb] * len(nz_y))
            comp_zi_col_list.extend([g_zi] * len(nz_y))

        # Phantom zero: one row with weight = count of zeros
        n_zeros = int(np.sum(~nz_mask))
        if n_zeros > 0:
            comp_y_list.append(np.array([0], dtype=np.int64))
            comp_w_list.append(np.array([float(n_zeros)]))
            comp_nb_col_list.append(g_nb)
            comp_zi_col_list.append(g_zi)

    comp_y = np.concatenate(comp_y_list).astype(np.int64)
    comp_w = np.concatenate(comp_w_list).astype(np.float64)
    n_c = len(comp_y)

    # Build sparse design matrices
    rows = np.arange(n_c)
    nb_cols_arr = np.array(comp_nb_col_list, dtype=np.int32)
    zi_cols_arr = np.array(comp_zi_col_list, dtype=np.int32)

    # NB: intercept (col 0) for all rows + treatment column for non-reference
    nb_row_idx = list(range(n_c))  # intercept entries
    nb_col_idx = [0] * n_c
    nb_data = [1.0] * n_c
    treatment_mask = nb_cols_arr > 0
    for r in np.where(treatment_mask)[0]:
        nb_row_idx.append(int(r))
        nb_col_idx.append(int(nb_cols_arr[r]))
        nb_data.append(1.0)
    comp_nb = sp.csr_matrix((nb_data, (nb_row_idx, nb_col_idx)), shape=(n_c, n_nb))

    # ZI: one-hot
    comp_zi = sp.csr_matrix(
        ([1.0] * n_c, (list(range(n_c)), [int(c) for c in zi_cols_arr])),
        shape=(n_c, n_zi)
    )

    assert abs(comp_w.sum() - n) < 0.5, f"Weight sum {comp_w.sum()} != {n}"
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
    rep_idx = np.argmax(Z_arr, axis=1)
    rep_labels = [zi_names[j].replace('C(rep_id)[', '').replace(']', '') for j in rep_idx]

    df = pd.DataFrame({'rep_id': rep_labels, 'zi': zi_pred})
    return df.groupby('rep_id')['zi'].mean()


# ── Main fitting loop ───────────────────────────────────────────────────────
phantom_results = {}
phantom_zi = {}
saved_zi = {}

total_t0 = time.time()

for ct in cell_types:
    design = all_designs[ct]
    saved = saved_ortho.model[ct]
    for wk in list(saved['weights'].keys()):
        saved['weights'][wk] = np.array(saved['weights'][wk])

    print(f"\n{'='*70}", flush=True)
    print(f"Cell type: {ct}", flush=True)
    print(f"  Full data: {design['regressand'].shape[0]:,} rows, "
          f"nb={design['nb_regressors'].shape[1]} cols", flush=True)

    t0 = time.time()
    comp_y, comp_w, comp_nb, comp_zi, n_full, n_comp = compress(design)
    t_compress = time.time() - t0
    print(f"  Compressed: {n_comp:,} rows ({n_full/n_comp:.0f}x) in {t_compress:.1f}s", flush=True)

    t0 = time.time()
    zinbo = TensorZINB(
        endog=comp_y, exog=comp_nb, exog_infl=comp_zi,
        sample_weight=comp_w, nb_only=False
    )
    res = zinbo.fit(return_history=True, init_method="nb", device_type="CPU")
    fit_time = time.time() - t0

    # Label weights
    nb_names = list(design['nb_regressors'].columns)
    zi_names = list(design['zi_regressors'].columns)
    res['weights']['x_mu'] = pd.Series(res['weights']['x_mu'].flatten(), index=nb_names)
    res['weights']['x_pi'] = pd.Series(res['weights']['x_pi'].flatten(), index=zi_names)

    phantom_results[ct] = res

    # ZI extraction
    phantom_zi[ct] = extract_zi(res['weights'], design)
    saved_zi[ct] = extract_zi(saved['weights'], design)

    # Compare
    llf_pct = 100 * (res['llf_total'] - saved['llf_total']) / abs(saved['llf_total'])

    # x_mu comparison
    s_mu = np.array(saved['weights']['x_mu']).flatten()
    p_mu = np.array(res['weights']['x_mu']).flatten()
    mu_pct = np.mean(np.abs((p_mu - s_mu) / np.where(np.abs(s_mu) > 1e-6, s_mu, 1.0)) * 100)

    # x_pi comparison
    s_pi = np.array(saved['weights']['x_pi']).flatten()
    p_pi = np.array(res['weights']['x_pi']).flatten()
    pi_pct = np.mean(np.abs((p_pi - s_pi) / np.where(np.abs(s_pi) > 1e-6, s_pi, 1.0)) * 100)

    print(f"  Fit: {fit_time:.1f}s, epochs={res['epochs']}", flush=True)
    print(f"  LLF %diff: {llf_pct:+.4f}%", flush=True)
    print(f"  x_mu mean |%diff|: {mu_pct:.2f}%", flush=True)
    print(f"  x_pi mean |%diff|: {pi_pct:.2f}%", flush=True)
    print(f"  ZI per rep (phantom): {dict(phantom_zi[ct].round(4))}", flush=True)
    print(f"  ZI per rep (saved):   {dict(saved_zi[ct].round(4))}", flush=True)

    del zinbo

total_time = time.time() - total_t0
print(f"\nAll {len(cell_types)} cell types fit in {total_time:.1f}s ({total_time/60:.1f}min)", flush=True)

# ── Save results ────────────────────────────────────────────────────────────
print("\nSaving results...", flush=True)
with open(f'{OUT_DIR}/phantom_by_cell_type_results.pkl', 'wb') as f:
    pickle.dump(phantom_results, f)

# ── Bounds-style ZI comparison ──────────────────────────────────────────────
print("\n" + "="*80, flush=True)
print("BOUNDS-STYLE by_cell_type ZI COMPARISON", flush=True)
print("="*80, flush=True)

# Average ZI across all cell types per rep
phantom_zi_dfs = [phantom_zi[ct].rename(ct) for ct in cell_types]
saved_zi_dfs = [saved_zi[ct].rename(ct) for ct in cell_types]

phantom_zi_all = pd.concat(phantom_zi_dfs, axis=1)
saved_zi_all = pd.concat(saved_zi_dfs, axis=1)

phantom_by_ct_zi = phantom_zi_all.mean(axis=1)
saved_by_ct_zi = saved_zi_all.mean(axis=1)

from scMPRAforge.core import SHENDURE_BOUNDS

print("\nPer-rep by_cell_type_zi comparison:", flush=True)
print(f"{'rep':<6} {'SHENDURE_BOUNDS':>16} {'Saved ortho':>12} {'Phantom':>12} {'%diff':>12}", flush=True)
print("-"*62, flush=True)

for rep in sorted(phantom_by_ct_zi.index):
    sb = SHENDURE_BOUNDS.by_cell_type_zi.get(rep, float('nan'))
    sv = saved_by_ct_zi.get(rep, float('nan'))
    ph = phantom_by_ct_zi.get(rep, float('nan'))
    pct = 100 * (ph - sv) / abs(sv) if abs(sv) > 1e-10 else 0
    print(f"{rep:<6} {sb:>16.6f} {sv:>12.6f} {ph:>12.6f} {pct:>+11.4f}%", flush=True)

print(f"\nMean across reps:", flush=True)
print(f"  SHENDURE_BOUNDS: {SHENDURE_BOUNDS.by_cell_type_zi.mean():.6f}", flush=True)
print(f"  Saved ortho:     {saved_by_ct_zi.mean():.6f}", flush=True)
print(f"  Phantom:         {phantom_by_ct_zi.mean():.6f}", flush=True)
pct_mean = 100*(phantom_by_ct_zi.mean() - saved_by_ct_zi.mean())/abs(saved_by_ct_zi.mean())
print(f"  %diff:           {pct_mean:+.4f}%", flush=True)

# Save comparison
comparison = {
    'phantom_by_ct_zi': phantom_by_ct_zi,
    'saved_by_ct_zi': saved_by_ct_zi,
    'shendure_bounds_by_ct_zi': SHENDURE_BOUNDS.by_cell_type_zi,
    'phantom_zi_all': phantom_zi_all,
    'saved_zi_all': saved_zi_all,
    'total_fit_time_s': total_time,
}
with open(f'{OUT_DIR}/bounds_comparison_by_cell_type.pkl', 'wb') as f:
    pickle.dump(comparison, f)

print(f"\nResults saved to {OUT_DIR}/", flush=True)
print("Done.", flush=True)
