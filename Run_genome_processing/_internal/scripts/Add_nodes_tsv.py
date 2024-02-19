import sys
import pandas as pd


def insert_nodes(reference_file, sorted_file, outliers_file, output_file):
    # Read reference TSV file into a DataFrame
    reference_df = pd.read_csv(reference_file, sep='\t')

    # Read sorted TSV file into a DataFrame
    sorted_df = pd.read_csv(sorted_file, sep='\t')

    # Read outliers TSV file into a DataFrame
    outliers_df = pd.read_csv(outliers_file, sep='\t')

    # Remove rows with no date values
    reference_df = reference_df.dropna(subset=['date'])
    sorted_df = sorted_df.dropna(subset=['date'])

    # Remove rows from sorted DataFrame where the first column matches outliers
    sorted_df = sorted_df[~sorted_df.iloc[:, 0].isin(outliers_df.iloc[:, 0])]

    # Insert NODE rows into sorted DataFrame based on numeric date column
    result_df = pd.concat([sorted_df, reference_df[reference_df['#node'].str.startswith('NODE_')]])

    # Removes rows with "--"
    result_df = result_df[result_df['date' or 'numeric date'] != '--']

    # Sort DataFrame by 'numeric date' column
    result_df.sort_values('numeric date', inplace=True)

    # Write resulting DataFrame to new TSV file
    result_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py reference_file.tsv sorted_file.tsv outliers_file.tsv output_file.tsv")
    else:
        reference_file = sys.argv[1]
        sorted_file = sys.argv[2]
        outliers_file = sys.argv[3]
        output_file = sys.argv[4]
        insert_nodes(reference_file, sorted_file, outliers_file, output_file)
