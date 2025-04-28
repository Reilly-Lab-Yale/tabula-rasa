import argparse
import pandas as pd


def generate_lookups(counts_file, output_file):
    # Load counts file
    df = pd.read_table(counts_file)

    # Check that each cre_id maps to exactly one cre_class
    mapping = df[['cre_id', 'cre_class']].drop_duplicates()

    # Find cre_ids that are associated with multiple cre_classes
    duplicated = mapping.groupby('cre_id').size()
    problematic = duplicated[duplicated > 1]

    if not problematic.empty:
        raise ValueError(f"Some cre_ids are associated with multiple cre_classes:\n{problematic}")

    # Otherwise, build the lookup strings
    lookups = mapping.apply(lambda row: f"{row['cre_id']}__{row['cre_class']}", axis=1)

    # Write to output file
    with open(output_file, 'w') as f:
        for lookup in lookups:
            f.write(f"{lookup}\n")
    print(f"Lookups saved to {output_file}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a list of lookup strings from counts file.")
    parser.add_argument('counts_file', type=str, help="Path to the counts file")
    parser.add_argument('output_file', type=str, help="Path to save the generated lookups")
    
    # Parse arguments
    args = parser.parse_args()

    # Call the function to generate the lookups
    generate_lookups(args.counts_file, args.output_file)

if __name__ == '__main__':
    main()