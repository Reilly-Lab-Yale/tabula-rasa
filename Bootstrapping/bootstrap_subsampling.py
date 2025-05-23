import argparse
import pandas as pd
import numpy as np
import os

def summarize_group(df, suffix=''):
    ## going forward from here mean_umis is always referring to gex normalized mpra umis 
    return pd.Series({
        f'n_integrations{suffix}': df[['cell_bc', 'transfection_bc']].drop_duplicates().shape[0],
        f'mean_umis_mpra_bc{suffix}': df['normalized_umis_mpra_bc'].mean()
    })

def bootstrap_sample(df_biol, idx_CRE_oi, idx_controls, n_integration_oi):
    sample_idx = np.concatenate(
        [np.random.choice(idxs, n_integration_oi, replace=True) for idxs in [idx_CRE_oi] + idx_controls]
    )
    return df_biol.loc[sample_idx]

def summarize_clusters(df_boot, group_cols, cell_types):
    # Add pseudocount to avoid zero division
    # df_boot['normalized_umis_mpra_bc'] = df_boot['normalized_umis_mpra_bc'] + 1

    # Now calculate expression in "not that cluster"
    summary = []

    for cell_type in cell_types:
        df_in_cluster = df_boot[df_boot['cell_type'] == cell_type]
        df_not_cluster = df_boot[df_boot['cell_type'] != cell_type]

        group_vals = df_in_cluster[group_cols].drop_duplicates()
        group_vals = group_vals.copy()
        group_vals['cell_type'] = cell_type  # preserve the cluster label

        temp_summary_not = df_not_cluster.groupby(['cre_id', 'cre_class', 'biol_rep']).apply(
            lambda g: summarize_group(g, suffix='_not_cluster'),
            include_groups=False
        ).reset_index()
        temp_summary = df_in_cluster.groupby(['cre_id', 'cre_class', 'biol_rep']).apply(
            lambda g: summarize_group(g, suffix=''),
            include_groups=False
        ).reset_index()
        temp_summary['cell_type'] = cell_type

        # Expand temp_summary to match cluster cell_type
        temp_summary = temp_summary_not.merge(temp_summary, 
                                            on=['cre_id', 'cre_class', 'biol_rep'], how='right')

        summary.append(temp_summary)

    summary_df = pd.concat(summary).reset_index(drop=True)

    # Compute fold change
    summary_df['FC_cluster_mBC_UMI'] = (
        summary_df['mean_umis_mpra_bc'] / summary_df['mean_umis_mpra_bc_not_cluster']
)

    return summary_df

def summarize_final(merged_summary):
    final_summary = merged_summary.groupby(['biol_rep', 'cre_id']).apply(
        lambda df: pd.Series({
            'max_expression_cluster_id': df.loc[df['mean_umis_mpra_bc'].idxmax(), 'cell_type'],
            'max_cluster_mBC_UMI': df['mean_umis_mpra_bc'].max(),
            'max_cluster_FC_mBC_UMI': df.loc[df['mean_umis_mpra_bc'].idxmax(), 'FC_cluster_mBC_UMI']
        }),
        include_groups=False
    ).reset_index()
    return final_summary

def load_biol_rep_mapping(mapping_path):
    mapping_df = pd.read_csv(mapping_path, sep='\s+', header=None, names=['rep_id', 'biol_rep'])
    print(mapping_df.head)
    return dict(zip(mapping_df['rep_id'], mapping_df['biol_rep']))

def main(lookup_str, n_bootstrap, date_str, counts_file, output_dir, control_cres, biol_rep_map_path):
    df = pd.read_table(counts_file)
    print("Columns loaded from counts file:", df.columns.tolist())

    if 'biol_rep' not in df.columns:
        print('mapping biol_rep')
        mapping_dict = load_biol_rep_mapping(biol_rep_map_path)
        df['biol_rep'] = df['rep_id'].map(mapping_dict)

    CRE_oi, CRE_type = lookup_str.split('__')

    biol_rep_list = sorted(df['biol_rep'].dropna().unique())
    cell_types = df['cell_type'].unique()

    print(f'biological replicate: {biol_rep_list}')
    group_cols = ['cre_id', 'cre_class', 'biol_rep', 'cell_type']

    # subset for CRE of interest + control CREs
    cre_set = [CRE_oi] + control_cres
    df_subset = df[df['cre_id'].isin(cre_set)]

    print(f"Running for CRE: {CRE_oi}")
    print(f"Subset size: {df_subset.shape[0]}")

    # Find integration numbers
    n_integrations = df_subset.groupby(['biol_rep', 'cre_id']).size().reset_index(name='n_int')
    median_integrations = df[df['cre_class'] == "devCRE"].groupby(['biol_rep', 'cre_id']).size().groupby('biol_rep').median().reset_index(name='median_int')
    median_integrations['median_int'] = median_integrations['median_int'].astype(int)

    results = []

    for biol_r in biol_rep_list:
        print(f'Processing biological replicate: {biol_r}')
        df_biol = df_subset[df_subset['biol_rep'] == biol_r]
        if df_biol.empty:
            print(f"No data for biological replicate {biol_r}, skipping.")
            continue

        idx_controls = [df_biol[df_biol['cre_id'] == ctrl].index for ctrl in control_cres]
        idx_CRE_oi = df_biol[df_biol['cre_id'] == CRE_oi].index
        if idx_CRE_oi.empty:
            print(f"No instances of {CRE_oi} in {biol_r}, skipping.")
            continue
        print(f"Biol rep: {biol_r}")
        print(f"  idx_CRE_oi size: {len(idx_CRE_oi)}")
        print(f"  idx_controls sizes: {[len(x) for x in idx_controls]}")

        if CRE_type == "devCRE":
            n_integration_oi = n_integrations.query('biol_rep == @biol_r and cre_id == @CRE_oi')['n_int'].values[0]
        else:
            n_integration_oi = median_integrations.query('biol_rep == @biol_r')['median_int'].values[0]
        print(f"  n_integration_oi: {n_integration_oi}")
        for boot in range(n_bootstrap):
            df_boot = bootstrap_sample(df_biol, idx_CRE_oi, idx_controls, n_integration_oi)
            merged_summary = summarize_clusters(df_boot, group_cols, cell_types) ## going forward from here mean_umis is always referring to gex normalized mpra umis 
            final_summary = summarize_final(merged_summary)

            final_summary['boot_id'] = boot + 1
            final_summary['bootstrap_for_CRE'] = CRE_oi
            final_summary['permuted'] = False
            results.append(final_summary)

    final_df = pd.concat(results)
    output_file = os.path.join(
        output_dir,
        f'bootstrap_subsampling_{CRE_type}_{CRE_oi}_{date_str}.txt'
    )
    final_df.to_csv(output_file, sep='\t', index=False)
    # Optional: check if any CREs are missing biological replicates
    if 'cre_id' in final_df.columns and 'biol_rep' in final_df.columns:
        cre_rep_counts = final_df.groupby('cre_id')['biol_rep'].nunique()
        expected_reps = final_df['biol_rep'].nunique()
        missing_reps = cre_rep_counts[cre_rep_counts < expected_reps]
        
        if not missing_reps.empty:
            print("Warning: Some CREs are missing biological replicates:")
            print(missing_reps)

    print(f'Saved final output to {output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Bootstrap subsampling for MPRA data.")
    parser.add_argument('--lookup_str', type=str, required=True, help="Lookup string in format CRE__type")
    parser.add_argument('--n_bootstrap', type=int, required=True, help="Number of bootstraps", default=10**4)
    parser.add_argument('--date_str', type=str, required=True, help="Date string for labeling outputs")
    parser.add_argument('--counts_file', type=str, required=True, help="Path to counts file (grouped counts)")
    parser.add_argument('--output_dir', type=str, required=True, help="Path to output folder")
    parser.add_argument('--control_cres', type=str, default="minP,noP", help="Comma-separated list of control CRE IDs")
    parser.add_argument('--biol_rep_map', type=str, required=True, help="Path to text file mapping rep_ids to biol_reps")

    args = parser.parse_args()

    main(
        lookup_str=args.lookup_str,
        n_bootstrap=args.n_bootstrap,
        date_str=args.date_str,
        counts_file=args.counts_file,
        output_dir=args.output_dir,
        control_cres=args.control_cres.split(','),
        biol_rep_map_path=args.biol_rep_map
    )