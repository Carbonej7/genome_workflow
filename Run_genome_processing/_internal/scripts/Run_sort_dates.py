import sys
from Bio import SeqIO
import pandas as pd

def main(fasta_file, metadata_file, sorted_metadata_output, reordered_sequences_output):
    # Phase 0: Check number of sequences
    header1 = []
    seq1 = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        header1.append(record.id)
        seq1.append(str(record.seq))

    # Read metadata
    T = pd.read_csv(metadata_file)

    M0 = len(T)
    M = len(seq1)
    N = len(seq1[0])

    if M0 != M:
        print("Warning: # of Sequence in CSV and FASTA files are not the same. Proceeding with available data...")

    # Phase 1: Sort the metadata
    try:
        T = T.sort_values(by='numeric date')
    except KeyError:
        print("Warning: 'numeric date column not found. Skipping sorting by date.")

    T.to_csv(sorted_metadata_output, index=False)  # Optionally save the sorted table back to CSV

    # Phase 2: Reordering sequence
    T = pd.read_csv(sorted_metadata_output)
    M0, _ = T.shape

    header2 = [''] * M0
    seq2 = [''] * M0
    flag2 = [True] * M0

    for k in range(M0):
        flag = 0

        for m in range(M):
            if flag2[k]:
                if T['#node'][k] == header1[m]:
                    header2[k] = header1[m]
                    seq2[k] = seq1[m]
                    flag = 1
                    flag2[k] = False
                    break

        if flag == 0:
            print(T['#node'][k], "is missing in the FASTA file. Skipping...")

    print("Number of sequences found:", sum(1 for x in header2 if x))

    # Write sequences to a new FASTA file
    with open(reordered_sequences_output, "w") as f:
        for i in range(M0):
            if header2[i] != '':
                f.write(">" + header2[i] + "\n")
                f.write(seq2[i] + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_fasta_file> <input_metadata_csv_file> <sorted_metadata_output> <reordered_sequences_output>")
        sys.exit(1)

    input_fasta_file = sys.argv[1]
    input_metadata_file = sys.argv[2]
    sorted_metadata_output = sys.argv[3]
    reordered_sequences_output = sys.argv[4]
    main(input_fasta_file, input_metadata_file, sorted_metadata_output, reordered_sequences_output)
