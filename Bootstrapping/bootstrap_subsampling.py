import argparse
import pandas as pd
import numpy as np
import os

def extract_true_uniques_from_str_series(series):
    all_values = []
    for val in series.dropna():
        if isinstance(val, str):
            parts = [v.strip() for v in val.split(',')]
            all_values.extend(parts)
    return np.unique(all_values)

def summarize_group(df, suffix=''):
    true_uniques = extract_true_uniques_from_str_series(df['transfection_bc'])
    return pd.Series({
        f'n_integrations{suffix}': len(true_uniques),
        f'mean_umis_mpra_bc{suffix}': df['umis_mpra_bc'].mean()
    })

def bootstrap_sample(df_biol, idx_CRE_oi, idx_controls, n_integration_oi):
    sample_idx = np.concatenate(
        [np.random.choice(idxs, n_integration_oi, replace=True) for idxs in [idx_CRE_oi] + idx_controls]
    )
    return df_biol.loc[sample_idx]

def summarize_clusters(df_boot, group_cols, cell_types):
    summary_cluster = df_boot.groupby(group_cols).apply(
        lambda g: summarize_group(g, suffix=''),
        include_groups=False
    ).reset_index()

    summary_not_cluster = []
    for cell_type in cell_types:
        df_not_cluster = df_boot[df_boot['cell_type'] != cell_type]
        temp_summary = df_not_cluster.groupby(group_cols).apply(
            lambda g: summarize_group(g, suffix='_not_cluster'),
            include_groups=False
        ).reset_index()
        temp_summary['cell_type'] = cell_type
        summary_not_cluster.append(temp_summary)

    summary_not_cluster_df = pd.concat(summary_not_cluster)
    merged_summary = summary_cluster.merge(summary_not_cluster_df, on=group_cols)
    merged_summary['FC_cluster_mBC_UMI'] = (
        merged_summary['mean_umis_mpra_bc'] / merged_summary['mean_umis_mpra_bc_not_cluster']
    )
    return merged_summary

def summarize_final(merged_summary):
    final_summary = merged_summary.groupby(['biol_rep', 'cre_id', 'cre_class']).apply(
        lambda df: pd.Series({
            'max_expression_cluster_id': df.loc[df['mean_umis_mpra_bc'].idxmax(), 'cell_type'],
            'max_cluster_mBC_UMI': df['mean_umis_mpra_bc'].max(),
            'max_cluster_FC_mBC_UMI': df.loc[df['mean_umis_mpra_bc'].idxmax(), 'FC_cluster_mBC_UMI']
        }),
        include_groups=False
    ).reset_index()
    return final_summary

def load_biol_rep_mapping(mapping_path):
    mapping_df = pd.read_csv(mapping_path, sep='\t', header=None, names=['rep_id', 'biol_rep'])
    return dict(zip(mapping_df['rep_id'], mapping_df['biol_rep']))

def main(lookup_str, n_bootstrap, date_str, counts_file, output_dir, control_cres, biol_rep_map_path):
    df = pd.read_table(counts_file)

    if 'biol_rep' not in df.columns:
        mapping_dict = load_biol_rep_mapping(biol_rep_map_path)
        df['biol_rep'] = df['rep_id'].map(mapping_dict)

    CRE_oi, CRE_type = lookup_str.split('__')

    biol_rep_list = sorted(df['biol_rep'].dropna().unique())
    cell_types = df['cell_type'].unique()

    group_cols = ['cre_id', 'cre_class', 'biol_rep', 'cell_type']

    # subset for CRE of interest + control CREs
    cre_set = [CRE_oi] + control_cres
    df_subset = df[df['cre_id'].isin(cre_set)]

    # Find integration numbers
    n_integrations = df_subset.groupby(['biol_rep', 'cre_id']).size().reset_index(name='n_int')
    median_integrations = df[df['cre_class'] == "devCRE"].groupby(['biol_rep', 'cre_id']).size().groupby('biol_rep').median().reset_index(name='median_int')
    median_integrations['median_int'] = median_integrations['median_int'].astype(int)

    results = []

    for biol_r in biol_rep_list:
        print(f'Processing biological replicate: {biol_r}')
        df_biol = df_subset[df_subset['biol_rep'] == biol_r]

        idx_CRE_oi = df_biol[df_biol['cre_id'] == CRE_oi].index
        idx_controls = [df_biol[df_biol['cre_id'] == ctrl].index for ctrl in control_cres]

        if CRE_type == "devCRE":
            n_integration_oi = n_integrations.query('biol_rep == @biol_r and cre_id == @CRE_oi')['n_int'].values[0]
        else:
            n_integration_oi = median_integrations.query('biol_rep == @biol_r')['median_int'].values[0]

        for boot in range(n_bootstrap):
            df_boot = bootstrap_sample(df_biol, idx_CRE_oi, idx_controls, n_integration_oi)
            merged_summary = summarize_clusters(df_boot, group_cols, cell_types)
            final_summary = summarize_final(merged_summary)

            final_summary['boot_id'] = boot + 1
            final_summary['bootstrap_for_CRE'] = CRE_oi
            results.append(final_summary)

    final_df = pd.concat(results)
    output_file = os.path.join(
        output_dir,
        f'bootstrap_subsampling_{CRE_type}_{CRE_oi}_{date_str}.txt'
    )
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f'Saved final output to {output_file}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Bootstrap subsampling for MPRA data.")
    parser.add_argument('lookup_str', type=str, help="Lookup string in format CRE__type__cluster_group")
    parser.add_argument('n_bootstrap', type=int, help="Number of bootstraps", default=10**4)
    parser.add_argument('date_str', type=str, help="Date string for labeling outputs")
    parser.add_argument('counts_file', type=str, help="Full path to counts file (e.g., shendure_counts_grouped.txt)")
    parser.add_argument('output_dir', type=str, help="Path to output folder")
    parser.add_argument('--control_cres', type=str, default="minP,noP", help="Comma-separated list of control CRE IDs")
    parser.add_argument('--biol_rep_map', type=str, required=False, help="Path to biol_rep map TSV file (rep_id<TAB>biol_rep)")

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