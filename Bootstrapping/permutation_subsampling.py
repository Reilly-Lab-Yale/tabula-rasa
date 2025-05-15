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

def permute_cell_types(df_biol):
    permuted_df = df_biol.copy()
    permuted_df['cell_type'] = np.random.permutation(permuted_df['cell_type'].values)
    return permuted_df

def summarize_clusters(df_perm, group_cols, cell_types):     ## going forward from here mean_umis is always referring to gex normalized mpra umis 

        # Add pseudocount to avoid zero division
    df_perm['normalized_umis_mpra_bc'] = df_perm['normalized_umis_mpra_bc'] + 1

    # Calculate expression in each cluster
    summary_cluster = df_perm[df_perm['cell_type'] == cell_type].groupby(group_cols).apply(
        lambda g: summarize_group(g, suffix=''),
        include_groups=False
    ).reset_index()

    # Now calculate expression in "not that cluster"
    summary_not_cluster = []
    for cell_type in cell_types:
        df_in_cluster = df_perm[df_perm['cell_type'] == cell_type]
        df_not_cluster = df_perm[df_perm['cell_type'] != cell_type]

        group_vals = df_in_cluster[group_cols].drop_duplicates()
        group_vals = group_vals.copy()
        group_vals['cell_type'] = cell_type  # preserve the cluster label

        temp_summary = df_not_cluster.groupby(['cre_id', 'cre_class', 'biol_rep']).apply(
            lambda g: summarize_group(g, suffix='_not_cluster'),
            include_groups=False
        ).reset_index()

        # Expand temp_summary to match cluster cell_type
        temp_summary = temp_summary.merge(group_vals[['cre_id', 'cre_class', 'biol_rep', 'cell_type']], 
                                          on=['cre_id', 'cre_class', 'biol_rep'], how='right')

        summary_not_cluster.append(temp_summary)

    summary_not_cluster_df = pd.concat(summary_not_cluster)

    # Merge cluster and not_cluster summaries
    merged_summary = summary_cluster.merge(summary_not_cluster_df, on=group_cols)

    # Compute fold change
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
    mapping_df = pd.read_csv(mapping_path, sep='\s+', header=None, names=['rep_id', 'biol_rep'])
    print(mapping_df.head())
    return dict(zip(mapping_df['rep_id'], mapping_df['biol_rep']))

def main(lookup_str, n_permutations, date_str, counts_file, output_dir, control_cres, biol_rep_map_path):
    df = pd.read_table(counts_file)
    print("Columns loaded from counts file:", df.columns.tolist())

    if 'biol_rep' not in df.columns:
        print('Mapping biol_rep')
        mapping_dict = load_biol_rep_mapping(biol_rep_map_path)
        df['biol_rep'] = df['rep_id'].map(mapping_dict)

    CRE_oi, CRE_type = lookup_str.split('__')

    biol_rep_list = sorted(df['biol_rep'].dropna().unique())
    cell_types = df['cell_type'].unique()

    print(f'Biological replicates: {biol_rep_list}')
    group_cols = ['cre_id', 'cre_class', 'biol_rep', 'cell_type']

    # Subset for CRE of interest + control CREs
    cre_set = [CRE_oi] + control_cres
    df_subset = df[df['cre_id'].isin(cre_set)]

    print(f"Running for CRE: {CRE_oi}")
    print(f"Subset size: {df_subset.shape[0]}")

    results = []
    for biol_r in biol_rep_list:
        print(f'Processing biological replicate: {biol_r}')
        df_biol = df_subset[df_subset['biol_rep'] == biol_r]
        if df_biol.empty:
            print(f"No data for biological replicate {biol_r}, skipping.")
            continue
        un_sum = summarize_clusters(df_biol, group_cols, cell_types)
        unperm_summary = summarize_final(un_sum)
        unperm_summary['perm_id'] = 'OBS'
        unperm_summary['permutation_for_CRE'] = CRE_oi
        results.append(unperm_summary)

        for perm in range(n_permutations):
            df_perm = permute_cell_types(df_biol)
            merged_summary = summarize_clusters(df_perm, group_cols, cell_types)     ## going forward from here mean_umis is always referring to gex normalized mpra umis 

            final_summary = summarize_final(merged_summary)

            final_summary['perm_id'] = perm + 1
            final_summary['permutation_for_CRE'] = CRE_oi
            results.append(final_summary)

    final_df = pd.concat(results)
    output_file = os.path.join(
        output_dir,
        f'permutation_subsampling_{CRE_type}_{CRE_oi}_{date_str}.txt'
    )
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f'Saved final output to {output_file}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Permutation-based subsampling for MPRA data.")
    parser.add_argument('--lookup_str', type=str, required=True, help="Lookup string in format CRE__type")
    parser.add_argument('--n_permutations', type=int, required=True, help="Number of permutations", default=10**4)
    parser.add_argument('--date_str', type=str, required=True, help="Date string for labeling outputs")
    parser.add_argument('--counts_file', type=str, required=True, help="Path to counts file (grouped counts)")
    parser.add_argument('--output_dir', type=str, required=True, help="Path to output folder")
    parser.add_argument('--control_cres', type=str, default="minP,noP", help="Comma-separated list of control CRE IDs")
    parser.add_argument('--biol_rep_map', type=str, required=True, help="Path to text file mapping rep_ids to biol_reps")

    args = parser.parse_args()

    main(
        lookup_str=args.lookup_str,
        n_permutations=args.n_permutations,
        date_str=args.date_str,
        counts_file=args.counts_file,
        output_dir=args.output_dir,
        control_cres=args.control_cres.split(','),
        biol_rep_map_path=args.biol_rep_map
    )