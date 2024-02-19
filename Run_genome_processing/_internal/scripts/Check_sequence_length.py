import sys
from Bio import SeqIO
import pandas as pd


def main(fasta_file, csv_file, output_fasta, output_csv):
    # Read the FASTA file
    fasta_records = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

    # Read the CSV file
    df = pd.read_csv(csv_file)

    # Check if the CSV file has a '#node' header
    if '#node' not in df.columns:
        print("CSV file doesn't have '#node' header.")
        exit()

    # Get the accession codes from the CSV file
    accession_codes = df['#node'].tolist()

    # Filter FASTA records based on accession codes
    filtered_records = [fasta_records[accession] for accession in accession_codes if accession in fasta_records]

    print("Number of sequences in original FASTA:", len(fasta_records))
    print("Number of sequences in original CSV:", len(accession_codes))

    # Write the filtered FASTA records to output file
    SeqIO.write(filtered_records, output_fasta, "fasta")

    # Update the CSV file to only contain accession codes present in the FASTA file
    removed_accessions = [accession for accession in accession_codes if accession not in fasta_records]
    df = df[df['#node'].isin(fasta_records)]
    df.to_csv(output_csv, index=False)

    # Print removed accessions
    if removed_accessions:
        print("Removed accessions:")
        for accession in removed_accessions:
            print(accession)

    # Print number of accessions in both outputs
    print("Number of sequences in filtered FASTA:", len(filtered_records))
    print("Number of sequences in filtered CSV:", len(df))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <fasta_file> <csv_file> <output_fasta> <output_csv>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    csv_file = sys.argv[2]
    output_fasta = sys.argv[3]
    output_csv = sys.argv[4]
    main(fasta_file, csv_file, output_fasta, output_csv)
