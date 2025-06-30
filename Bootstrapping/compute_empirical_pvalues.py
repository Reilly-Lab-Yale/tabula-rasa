#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy.stats import rankdata
import argparse
from statsmodels.stats.multitest import multipletests
from collections import Counter

def calculate_empirical_pvals(observed, null):
    """
    Calculate empirical p-value: proportion of null >= observed (plus 1).
    """
    pvals = np.array([(np.sum(null >= obs) + 1) / (len(null) + 1) for obs in observed])
    return pvals

def correct_pvals(pvals, method="fdr_bh"):
    """
    Apply multiple testing correction.
    Default is Benjamini-Hochberg (FDR).
    """
    corrected_pvals = multipletests(pvals, method=method)[1]
    return corrected_pvals

def main():
    parser = argparse.ArgumentParser(description="Calculate empirical p-values and FDR for activity and specificity")
    parser.add_argument("--bootstrap_file", required=True, help="Path to combined bootstrap results")
    parser.add_argument("--permutation_file", required=True, help="Path to combined permutation results")
    parser.add_argument("--output_activity_file", required=True, help="Path to output activity p-value table")
    parser.add_argument("--output_specificity_file", required=True, help="Path to output specificity p-value table")
    args = parser.parse_args()

    # Load
    df_boot = pd.read_csv(args.bootstrap_file, sep="\t")
    df_perm = pd.read_csv(args.permutation_file, sep="\t")

    # Determine CRE list
    CRE_list = df_boot['cre_id'].dropna().unique()
    biol_reps = df_boot['biol_rep'].dropna().unique()

    ### === ACTIVITY P-VALUES === ###
    activity_results = []

    for cre in CRE_list:
        for rep in biol_reps:
            df_boot_oi = df_boot[(df_boot['bootstrap_for_CRE'] == cre) & (df_boot['biol_rep'] == rep)]

            if df_boot_oi.empty:
                print(f"Warning: Missing data for cre {cre}, rep {rep}")
                continue
                   
            # expression values
            exp_cre = df_boot_oi[df_boot_oi['cre_id'] == cre]['max_cluster_mBC_UMI']
            exp_min_noP = df_boot_oi[df_boot_oi['bootstrap_for_CRE'] == cre][df_boot_oi['cre_id'] != cre]['max_cluster_mBC_UMI']


            if exp_cre.empty or exp_min_noP.empty:
                continue

            # empirical p-value
            empirical_pvals = calculate_empirical_pvals(exp_cre.values, exp_min_noP.values)

            # Most frequent cluster
            max_clusters = df_boot_oi[df_boot_oi['cre_id'] == cre]['max_expression_cluster_id']
            if not max_clusters.empty:
                most_common_cluster = Counter(max_clusters).most_common(1)[0][0]
                freq_most_common = np.mean(max_clusters == most_common_cluster)
            else:
                most_common_cluster = np.nan
                freq_most_common = np.nan

            activity_results.append({
                "cre_id": cre,
                "biol_rep": rep,
                "q25_bootstrap_max_exp_norm": np.quantile(exp_cre, 0.25),
                "q50_bootstrap_max_exp_norm": np.quantile(exp_cre, 0.5),
                "q75_bootstrap_max_exp_norm": np.quantile(exp_cre, 0.75),
                "empirical_bootstrap_p_val_vs_noP_minP": np.mean(empirical_pvals),
                "most_frequent_max_exp_cluster": most_common_cluster,
                "frequency_most_frequent_max_exp_cluster": freq_most_common
            })

    df_activity = pd.DataFrame(activity_results)

    # FDR correction by biological replicate
    df_activity["empirical_bootstrap_p_val_vs_noP_minP_BH_corr"] = np.nan
    for rep in biol_reps:
        mask = df_activity["biol_rep"] == rep
        if mask.sum() > 0:
            corrected = correct_pvals(df_activity.loc[mask, "empirical_bootstrap_p_val_vs_noP_minP"])
            df_activity.loc[mask, "empirical_bootstrap_p_val_vs_noP_minP_BH_corr"] = corrected

    # Save
    df_activity.to_csv(args.output_activity_file, sep="\t", index=False)

    print(f"Saved activity p-values to {args.output_activity_file}")

    # ### === CELL TYPE SPECIFICITY P-VALUES === ###
    # specificity_results = []

    # for cre in CRE_list:
    #     for rep in biol_reps:
    #         df_boot_oi = df_boot[(df_boot['bootstrap_for_CRE'] == cre) & (df_boot['biol_rep'] == rep) & (df_boot['cre_id'] == cre)]
    #         df_perm_oi = df_perm[(df_perm['permutation_for_CRE'] == cre) & (df_perm['biol_rep'] == rep) & (df_perm['cre_id'] == cre) & (df_perm['perm_id'] != 'OBS')]

    #         if df_perm_oi.empty:
    #             print(f"Warning: Missing data for cre {cre}, rep {rep}")
    #             continue

    #         # FC (fold change) values
    #         fc_boot = df_boot_oi['max_cluster_FC_mBC_UMI'].fillna(1)
    #         fc_perm = df_perm_oi['max_cluster_FC_mBC_UMI'].fillna(1)

    #         # empirical p-value
    #         empirical_pvals = calculate_empirical_pvals(fc_boot.values, fc_perm.values)

    #         # Most frequent cluster
    #         max_clusters = df_boot_oi['max_expression_cluster_id']
    #         if not max_clusters.empty:
    #             most_common_cluster = Counter(max_clusters).most_common(1)[0][0]
    #             freq_most_common = np.mean(max_clusters == most_common_cluster)
    #         else:
    #             most_common_cluster = np.nan
    #             freq_most_common = np.nan

    #         specificity_results.append({
    #             "cre_id": cre,
    #             "biol_rep": rep,
    #             "q25_bootstrap_FC_norm": np.quantile(fc_boot, 0.25),
    #             "q50_bootstrap_FC_norm": np.quantile(fc_boot, 0.5),
    #             "q75_bootstrap_FC_norm": np.quantile(fc_boot, 0.75),
    #             "empirical_permute_gex_cluster_p_val": np.mean(empirical_pvals),
    #             "most_frequent_max_exp_cluster": most_common_cluster,
    #             "frequency_most_frequent_max_exp_cluster": freq_most_common
    #         })

    # df_specificity = pd.DataFrame(specificity_results)

    # # FDR correction by biological replicate
    # df_specificity["empirical_permute_gex_cluster_p_val_BH_corr"] = np.nan
    # for rep in biol_reps:
    #     mask = df_specificity["biol_rep"] == rep
    #     if mask.sum() > 0:
    #         corrected = correct_pvals(df_specificity.loc[mask, "empirical_permute_gex_cluster_p_val"])
    #         df_specificity.loc[mask, "empirical_permute_gex_cluster_p_val_BH_corr"] = corrected

    # # Save
    # df_specificity.to_csv(args.output_specificity_file, sep="\t", index=False)

    # print(f"Saved specificity p-values to {args.output_specificity_file}")

if __name__ == "__main__":
    main()