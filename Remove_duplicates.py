import sys
from Bio import SeqIO
import pandas as pd


# This Code will read both fasta and metadata and check for duplicates and remove them if there are any
def check_duplicates(file_path):
    seen = set()
    duplicates = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() in seen:
                duplicates.append(line.strip())
            else:
                seen.add(line.strip())
    return duplicates


def filter_fasta_by_metadata(fasta_file, metadata_file, fasta_output_file, metadata_output_file):
    try:
        # Check for duplicates in metadata file
        metadata_duplicates = check_duplicates(metadata_file)
        if metadata_duplicates:
            print("Warning: Duplicates found in metadata file. Removing duplicates.")
            metadata_df = pd.read_csv(metadata_file, sep='\t')
            metadata_df.drop_duplicates(inplace=True)
            metadata_df.to_csv(metadata_output_file, sep='\t', index=False)
        else:
            metadata_df = pd.read_csv(metadata_file, sep='\t')
            metadata_df.to_csv(metadata_output_file, sep='\t', index=False)

        # Check for duplicates in FASTA file
        fasta_duplicates = check_duplicates(fasta_file)
        if fasta_duplicates:
            print("Warning: Duplicates found in FASTA file. Removing duplicates.")
            unique_records = {}
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if record.seq not in unique_records:
                    unique_records[record.seq] = record
            SeqIO.write(unique_records.values(), fasta_output_file, 'fasta')
        else:
            SeqIO.write(SeqIO.parse(fasta_file, 'fasta'), fasta_output_file, 'fasta')

        print("Filtered files created successfully.")

    except FileNotFoundError:
        print("Error: One of the input files not found.")
    except pd.errors.ParserError:
        print("Error: Unable to parse the metadata file.")
    except Exception as e:
        print("An error occurred:", str(e))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python Remove_duplicates.py fasta_file metadata_file fasta_output_file metadata_output_file")
        sys.exit(1)

    fasta_file = sys.argv[1]
    metadata_file = sys.argv[2]
    fasta_output_file = sys.argv[3]
    metadata_output_file = sys.argv[4]

    filter_fasta_by_metadata(fasta_file, metadata_file, fasta_output_file, metadata_output_file)
