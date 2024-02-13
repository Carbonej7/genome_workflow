import sys
import re
import pandas as pd


def clean_metadata(input_metadata_file, cleaned_metadata_file):
    # Read the metadata CSV file into a pandas DataFrame
    metadata_df = pd.read_csv(input_metadata_file)

    # Check if there are any columns other than "Accession Code" and "Collection_Date"
    other_columns = set(metadata_df.columns) - {'Accession', 'Collection_Date'}
    if other_columns:
        raise ValueError(
            f"Found unexpected columns: {', '.join(other_columns)}. Only 'Accession' and 'Collection_Date' are"
            f"allowed.")

    # Keep only the "Accession Code" and "Collection_Date" columns
    metadata_df = metadata_df[['Accession', 'Collection_Date']]

    # Keep track of encountered accession codes
    encountered_accession_codes = set()

    # Iterate through each row in the DataFrame
    cleaned_metadata_rows = []
    for index, row in metadata_df.iterrows():
        accession_code = row['Accession']
        # If the accession code is not a duplicate, add it to the cleaned metadata
        if accession_code not in encountered_accession_codes:
            cleaned_metadata_rows.append(row)
            encountered_accession_codes.add(accession_code)

    # Create a new DataFrame with the cleaned metadata rows
    cleaned_metadata_df = pd.DataFrame(cleaned_metadata_rows)

    # Write the cleaned metadata to a new CSV file
    cleaned_metadata_df.to_csv(cleaned_metadata_file, index=False)

    return cleaned_metadata_file


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python Run_clean_metadata.py input_metadata.csv output_clean_metadata.csv")
        sys.exit(1)

    input_metadata_file = sys.argv[1]
    cleaned_metadata_file = sys.argv[2]
    clean_metadata(input_metadata_file, cleaned_metadata_file)
    print(f'Cleaned version numbers from metadata saved to: {cleaned_metadata_file}')
