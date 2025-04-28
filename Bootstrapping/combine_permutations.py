import argparse
import pandas as pd
import glob
import os

def main(input_dir, output_file):
    # Find all permutation output files
    pattern = os.path.join(input_dir, 'permutation_subsampling_*_*.txt')
    files = glob.glob(pattern)
    print(f"Found {len(files)} files to combine.")

    if not files:
        raise ValueError(f"No permutation files found in {input_dir}")

    # Read and concatenate all permutation DataFrames
    dfs = []
    for file in files:
        try:
            df = pd.read_csv(file, sep='\t')
            dfs.append(df)
        except Exception as e:
            print(f"Failed to read {file}: {e}")

    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Combined dataframe shape: {combined_df.shape}")

    # Save combined dataframe
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Combined file saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine permutation output files into one file.")
    parser.add_argument('--input_dir', type=str, required=True, help='Input directory containing permutation files')
    parser.add_argument('--output_file', type=str, required=True, help='Path to save combined output')
    args = parser.parse_args()

    main(args.input_dir, args.output_file)