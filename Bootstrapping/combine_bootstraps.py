import argparse
import pandas as pd
import os
import glob

def combine_all_bootstrap_tables(input_dir, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Looking for bootstrap_subsampling files in {input_dir}...")

    # Only grab files matching the bootstrap_subsampling_*_*.txt pattern
    pattern = os.path.join(input_dir, "bootstrap_subsampling_*_*.txt")
    files = glob.glob(pattern)

    if not files:
        print("No matching bootstrap files found!")
        return

    print(f"Found {len(files)} bootstrap files to combine.")

    # Read first file to get expected columns
    df_example = pd.read_csv(files[0], sep="\t")
    n_rows_expected = df_example.shape[0]
    columns = df_example.columns

    all_dfs = []
    problem_files = []

    for file in files:
        df = pd.read_csv(file, sep="\t")
        if df.shape[1] != len(columns):
            print(f"Warning: {file} has different number of columns (expected {len(columns)})")
            problem_files.append(file)
            continue
        all_dfs.append(df)

    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        combined_df.to_csv(output_file, sep="\t", index=False)
        print(f"Saved combined file to {output_file}")
    else:
        print("No valid bootstrap files to combine.")

    if problem_files:
        print("Files with unexpected number of columns:")
        for pf in problem_files:
            print(f"  - {pf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine bootstrap output tables.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing bootstrap output files.")
    parser.add_argument('--output_file', type=str, required=True, help="Path for the combined output file.")
    args = parser.parse_args()

    combine_all_bootstrap_tables(args.input_dir, args.output_file)