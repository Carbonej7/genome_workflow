import sys
from Bio import SeqIO
import pandas as pd


def read_fasta(file_path):
    headers = []
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
    print("Read", len(headers), "sequences from the FASTA file")
    return headers, sequences


def read_tsv(file_path):
    return pd.read_csv(file_path, sep='\t')


def write_tsv(data, file_path):
    data.to_csv(file_path, sep='\t', index=False)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "Usage: python Run_sort_dates_2.py metadata_file.tsv fasta_file.fasta output_metadata.tsv "
            "output_fasta.fasta")
        sys.exit(1)

    metadata_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_metadata = sys.argv[3]
    output_fasta = sys.argv[4]

    # Phase 0: Check number of sequences
    headers, sequences = read_fasta(fasta_file)
    print("Number of sequences in the FASTA file:", len(headers))

    if not sequences:
        print("Error: No sequences found in the provided FASTA file.")
        sys.exit(1)

    # Read and sort the TSV file based on the numeric date column
    try:
        T = read_tsv(metadata_file)
        T = T.sort_values(by='numeric date')
        write_tsv(T, output_metadata)
    except KeyError:
        print("Error: 'numeric date' column not found in the TSV file.")
        sys.exit(1)

    # Reorder sequences based on the sorted TSV file
    ordered_headers = list(T['#node'])
    ordered_sequences = []
    for header in ordered_headers:
        if header in headers:
            index = headers.index(header)
            ordered_sequences.append(sequences[index])
        else:
            print(f"Warning: Header '{header}' from TSV file not found in FASTA file. Skipping...")
            # You can choose to handle missing headers as per your requirement

    # Write sequences to a new FASTA file
    with open(output_fasta, "w") as f:
        for header, sequence in zip(ordered_headers, ordered_sequences):
            f.write(">" + header + "\n")
            f.write(sequence + "\n")
