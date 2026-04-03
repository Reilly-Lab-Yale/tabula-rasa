"""
Test: Shendure by_cell_type CM fit using phantom-zero compression path on GPU.

Fits by_cell_type models through the standard_fit CM path (which now uses
_build_cm_fit_inputs for phantom-zero compression) and compares to the saved
CM ZINB ortho on disk.

Note: _mom_from_cm_maps is a stub that returns None, so these fits use NB init
instead of MoM.  The saved ortho used MoM init, so per-model parameters will
differ more than the by_cre test.  Bounds-level ZI should still be comparable.
"""
import pickle, sys, time, os, subprocess, json, re
import numpy as np
import pandas as pd

sys.path.insert(0, '/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa-cohen-regen')

from pathlib import Path

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data/shendure")
SAVED_ORTHO = DATA_DIR / "shendure_ortho_consider_missing_20260320"
SLURM_JOB_ID = os.environ.get("SLURM_JOB_ID", "unknown")


def main():
    from dask.distributed import Client, LocalCluster, performance_report
    import scMPRAforge.core as scm

    # ── Setup Dask (GPU-capable) ───────────────────────────────────────
    print("Setting up Dask local cluster with GPU resource...", flush=True)
    cluster = LocalCluster(
        n_workers=1,
        threads_per_worker=1,
        processes=False,
        memory_limit="60GB",
    )
    client = Client(cluster)
    print(f"Dashboard: {client.dashboard_link}", flush=True)

    # ── Load data ──────────────────────────────────────────────────────
    print("Loading Shendure data...", flush=True)
    shendure = scm.scMPRA_data.from_tsv(str(DATA_DIR / "shendure_processed.tsv"))
    shendure.set_negative_controls(["minP", "noP"])
    shendure.set_reference_cell("Pluripotent")
    shendure.ortho_filter()
    shendure.set_consider_missing(True)
    print("Data loaded and filtered.", flush=True)

    # ── Load saved ortho for comparison ────────────────────────────────
    print(f"\nLoading saved CM ortho from {SAVED_ORTHO}...", flush=True)
    with open(SAVED_ORTHO / "by_cell_type.pkl", "rb") as f:
        saved_by_ct = pickle.load(f)

    saved_cts = sorted(saved_by_ct.model.keys())
    print(f"Saved ortho has {len(saved_cts)} cell types", flush=True)

    for ct in saved_cts:
        for wk in list(saved_by_ct.model[ct]['weights'].keys()):
            saved_by_ct.model[ct]['weights'][wk] = np.array(
                saved_by_ct.model[ct]['weights'][wk])

    # ── Fit by_cell_type with phantom-zero CM + GPU ────────────────────
    print("\nFitting by_cell_type with phantom-zero CM path (GPU)...", flush=True)
    ortho_obj = scm.ortho()

    report_path = str(SCRIPT_DIR / "dask_performance_report.html")
    t0 = time.time()
    with performance_report(filename=report_path):
        ortho_obj.fit_by_cell_type_models(
            client=client,
            dat=shendure,
            fit_resources={},
        )
    fit_time = time.time() - t0
    print(f"\nFit complete in {fit_time:.1f}s ({fit_time/60:.1f}min)", flush=True)
    print(f"Dask performance report: {report_path}", flush=True)

    # ── Gather results ─────────────────────────────────────────────────
    print("Gathering results...", flush=True)
    by_ct = ortho_obj.by_cell_type
    for ct in list(by_ct.model.keys()):
        fut = by_ct.model[ct]
        by_ct.model[ct] = fut.result() if hasattr(fut, 'result') else fut

    phantom_cts = sorted(by_ct.model.keys())
    print(f"Phantom fit has {len(phantom_cts)} cell types", flush=True)

    # ── Compare ────────────────────────────────────────────────────────
    print("\n" + "=" * 80, flush=True)
    print("COMPARISON: phantom-zero CM vs saved CM ortho (by_cell_type)", flush=True)
    print("=" * 80, flush=True)

    common_cts = sorted(set(saved_cts) & set(phantom_cts))
    print(f"Comparing {len(common_cts)} common cell types "
          f"(saved={len(saved_cts)}, phantom={len(phantom_cts)})", flush=True)
    print("NOTE: saved ortho used MoM init; phantom uses NB init (stub).",
          flush=True)

    llf_diffs = []
    mu_diffs = []
    pi_diffs = []
    theta_diffs = []

    for ct in common_cts:
        saved = saved_by_ct.model[ct]
        phantom = by_ct.model[ct]

        for wk in list(phantom['weights'].keys()):
            phantom['weights'][wk] = np.array(phantom['weights'][wk])

        llf_pct = 100 * (phantom['llf_total'] - saved['llf_total']) / abs(saved['llf_total'])
        llf_diffs.append(llf_pct)

        s_mu = np.array(saved['weights']['x_mu']).flatten()
        p_mu = np.array(phantom['weights']['x_mu']).flatten()
        denom = np.where(np.abs(s_mu) > 1e-6, s_mu, 1.0)
        mu_pct = np.mean(np.abs((p_mu - s_mu) / denom) * 100)
        mu_diffs.append(mu_pct)

        if 'x_pi' in saved['weights'] and 'x_pi' in phantom['weights']:
            s_pi = np.array(saved['weights']['x_pi']).flatten()
            p_pi = np.array(phantom['weights']['x_pi']).flatten()
            denom = np.where(np.abs(s_pi) > 1e-6, s_pi, 1.0)
            pi_pct = np.mean(np.abs((p_pi - s_pi) / denom) * 100)
            pi_diffs.append(pi_pct)

        s_th = float(np.array(saved['weights']['theta']).flatten()[0])
        p_th = float(np.array(phantom['weights']['theta']).flatten()[0])
        th_pct = 100 * abs(p_th - s_th) / abs(s_th) if abs(s_th) > 1e-10 else 0
        theta_diffs.append(th_pct)

    llf_diffs = np.array(llf_diffs)
    mu_diffs = np.array(mu_diffs)
    pi_diffs = np.array(pi_diffs)
    theta_diffs = np.array(theta_diffs)

    print(f"\nLLF %diff:   mean={llf_diffs.mean():+.4f}%, "
          f"std={llf_diffs.std():.4f}%, max |%diff|={np.abs(llf_diffs).max():.4f}%",
          flush=True)
    print(f"x_mu |%diff|: mean={mu_diffs.mean():.2f}%, "
          f"median={np.median(mu_diffs):.2f}%, max={mu_diffs.max():.2f}%",
          flush=True)
    if len(pi_diffs) > 0:
        print(f"x_pi |%diff|: mean={pi_diffs.mean():.2f}%, "
              f"median={np.median(pi_diffs):.2f}%, max={pi_diffs.max():.2f}%",
              flush=True)
    print(f"theta |%diff|: mean={theta_diffs.mean():.2f}%, "
          f"median={np.median(theta_diffs):.2f}%, max={theta_diffs.max():.2f}%",
          flush=True)

    # Per-cell-type detail
    print(f"\n{'Cell type':<20} {'Saved LLF':>14} {'Phantom LLF':>14} {'LLF %diff':>12}",
          flush=True)
    print("-" * 64, flush=True)
    for ct in common_cts:
        sl = saved_by_ct.model[ct]['llf_total']
        pl = by_ct.model[ct]['llf_total']
        pct = 100 * (pl - sl) / abs(sl)
        print(f"{ct:<20} {sl:>14.2f} {pl:>14.2f} {pct:>+11.4f}%", flush=True)

    # ── ZI comparison (bounds-style) ──────────────────────────────────
    print("\n" + "=" * 80, flush=True)
    print("BOUNDS-STYLE ZI COMPARISON", flush=True)
    print("=" * 80, flush=True)

    print("Loading saved design matrices for ZI extraction...", flush=True)
    saved_design_dir = SAVED_ORTHO / "by_cell_type_design"
    with open(saved_design_dir / "_keys.json") as f:
        design_keys = json.load(f)
    all_designs = {}
    for i, k in enumerate(design_keys):
        with open(saved_design_dir / f"{i}.pkl", "rb") as f:
            all_designs[k] = pickle.load(f)

    def extract_zi(model_weights, design_mats):
        import scipy.sparse as sp
        Z = design_mats['zi_regressors']
        zi_names = design_mats.get('zi_regressor_names',
                                   list(Z.columns) if hasattr(Z, 'columns') else [])
        Z_arr = Z.toarray() if sp.issparse(Z) else np.array(Z)
        x_pi = np.array(model_weights['x_pi']).flatten()
        linear_zi = Z_arr @ x_pi
        zi_pred = 1 / (1 + np.exp(-linear_zi))
        rep_idx = np.argmax(Z_arr, axis=1)
        rep_labels = []
        for j in rep_idx:
            name = zi_names[j] if j < len(zi_names) else str(j)
            if '[' in name:
                m = re.search(r'\[(.*?)\]', name)
                name = m.group(1) if m else name
            rep_labels.append(name)
        df = pd.DataFrame({'rep_id': rep_labels, 'zi': zi_pred})
        return df.groupby('rep_id')['zi'].mean()

    phantom_zi_per_ct = {}
    saved_zi_per_ct = {}

    for ct in common_cts:
        if ct not in all_designs:
            continue
        design = all_designs[ct]
        if 'x_pi' in by_ct.model[ct]['weights']:
            phantom_zi_per_ct[ct] = extract_zi(by_ct.model[ct]['weights'], design)
        if 'x_pi' in saved_by_ct.model[ct]['weights']:
            saved_zi_per_ct[ct] = extract_zi(saved_by_ct.model[ct]['weights'], design)

    if phantom_zi_per_ct and saved_zi_per_ct:
        phantom_zi_all = pd.concat(
            [phantom_zi_per_ct[c].rename(c) for c in phantom_zi_per_ct], axis=1)
        saved_zi_all = pd.concat(
            [saved_zi_per_ct[c].rename(c) for c in saved_zi_per_ct], axis=1)

        phantom_by_ct_zi = phantom_zi_all.mean(axis=1)
        saved_by_ct_zi = saved_zi_all.mean(axis=1)

        print("\nPer-rep by_cell_type_zi comparison:", flush=True)
        print(f"{'rep':<6} {'Saved CM':>12} {'Phantom CM':>12} {'%diff':>12}", flush=True)
        print("-" * 46, flush=True)
        for rep in sorted(phantom_by_ct_zi.index):
            sv = saved_by_ct_zi.get(rep, float('nan'))
            ph = phantom_by_ct_zi.get(rep, float('nan'))
            pct = 100 * (ph - sv) / abs(sv) if abs(sv) > 1e-10 else 0
            print(f"{rep:<6} {sv:>12.6f} {ph:>12.6f} {pct:>+11.4f}%", flush=True)

        print(f"\nMean across reps:", flush=True)
        print(f"  Saved CM:    {saved_by_ct_zi.mean():.6f}", flush=True)
        print(f"  Phantom CM:  {phantom_by_ct_zi.mean():.6f}", flush=True)
        pct_mean = 100 * (phantom_by_ct_zi.mean() - saved_by_ct_zi.mean()) / abs(saved_by_ct_zi.mean())
        print(f"  %diff:       {pct_mean:+.4f}%", flush=True)

        common_zi_cts = sorted(set(phantom_zi_per_ct.keys()) & set(saved_zi_per_ct.keys()))
        phantom_ct_means = phantom_zi_all.mean(axis=0)
        saved_ct_means = saved_zi_all.mean(axis=0)
        common_cols = phantom_ct_means.index.intersection(saved_ct_means.index)
        if len(common_cols) > 2:
            r = np.corrcoef(phantom_ct_means[common_cols],
                            saved_ct_means[common_cols])[0, 1]
            print(f"\nPer-cell-type ZI Pearson r: {r:.4f}", flush=True)
        print(f"Per-cell-type ZI std: saved={saved_ct_means[common_cols].std():.6f}, "
              f"phantom={phantom_ct_means[common_cols].std():.6f}", flush=True)

    # ── Save comparison pickle ─────────────────────────────────────────
    comparison = {
        'llf_diffs': llf_diffs,
        'mu_diffs': mu_diffs,
        'pi_diffs': pi_diffs,
        'theta_diffs': theta_diffs,
        'fit_time_s': fit_time,
        'n_cell_types': len(common_cts),
    }
    if phantom_zi_per_ct and saved_zi_per_ct:
        comparison['phantom_by_ct_zi'] = phantom_by_ct_zi
        comparison['saved_by_ct_zi'] = saved_by_ct_zi

    with open(SCRIPT_DIR / 'comparison_results_by_cell_type.pkl', 'wb') as f:
        pickle.dump(comparison, f)

    print(f"\nResults saved to {SCRIPT_DIR}/comparison_results_by_cell_type.pkl",
          flush=True)

    client.close()
    cluster.close()

    # ── Collect resource usage stats ───────────────────────────────────
    print("\n" + "=" * 80, flush=True)
    print("RESOURCE USAGE STATISTICS", flush=True)
    print("=" * 80, flush=True)

    if SLURM_JOB_ID != "unknown":
        sacct_fields = (
            "JobID,JobName,State,Partition,NodeList,AllocCPUS,ReqMem,"
            "Submit,Start,End,Elapsed,TotalCPU,UserCPU,SystemCPU,"
            "MaxRSS,AveRSS,MaxVMSize,ReqTRES,AllocTRES"
        )
        sacct_out = subprocess.run(
            ["sacct", "-j", SLURM_JOB_ID, f"--format={sacct_fields}", "--parsable2"],
            capture_output=True, text=True
        ).stdout
        print(f"\nsacct for job {SLURM_JOB_ID}:", flush=True)
        print(sacct_out, flush=True)

        # Get node hardware info
        node_name = os.environ.get("SLURMD_NODENAME", "")
        if node_name:
            scontrol_out = subprocess.run(
                ["scontrol", "show", "node", node_name],
                capture_output=True, text=True
            ).stdout
            # Extract key fields
            for field in ["CPUTot", "RealMemory", "AvailableFeatures",
                          "Gres", "CfgTRES", "NodeName", "Arch"]:
                for line in scontrol_out.split():
                    if line.startswith(f"{field}="):
                        print(f"  {line}", flush=True)

            # GPU info via nvidia-smi
            gpu_info = subprocess.run(
                ["nvidia-smi", "--query-gpu=name,memory.total,driver_version,compute_cap",
                 "--format=csv,noheader"],
                capture_output=True, text=True
            ).stdout.strip()
            if gpu_info:
                print(f"\nGPU: {gpu_info}", flush=True)

        # jobstats if available
        jobstats_out = subprocess.run(
            ["jobstats", "--no-color", SLURM_JOB_ID],
            capture_output=True, text=True
        )
        if jobstats_out.returncode == 0 and jobstats_out.stdout.strip():
            print(f"\njobstats:\n{jobstats_out.stdout}", flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
