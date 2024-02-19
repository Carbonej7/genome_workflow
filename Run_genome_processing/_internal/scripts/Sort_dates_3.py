import pandas as pd
import sys
from Bio import SeqIO
import tempfile


def add_sequences(reference_file, input_file):
    # Read fasta and read accessions "NODE"
    reference_sequences = {}
    for record in SeqIO.parse(reference_file, "fasta"):
        if record.id.startswith("NODE"):
            reference_sequences[record.id] = record.seq

    # Read input
    input_sequences = []
    for record in SeqIO.parse(input_file, "fasta"):
        input_sequences.append(record)

    # Add reference "NODE" sequences to input
    for accession, sequence in reference_sequences.items():
        input_sequences.append(SeqIO.SeqRecord(seq=sequence, id=accession, description=""))
    # Creation of temporary output
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_sequences:
        SeqIO.write(input_sequences, temp_sequences, "fasta")
        temp_file_name = temp_sequences.name

    return temp_file_name


def sort_sequences_by_date(fasta_file, csv_file, output_file):
    # Use SeqIO to parse list of fasta - sequences not dictionary, so must convert list
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    csv_data = pd.read_csv(csv_file)
    # Make a dictionary to map columns #node and "numeric date
    node_to_date = dict(zip(csv_data['#node'], csv_data['numeric date']))
    # Sort sequences by numeric date in dictionary, float('inf') sorts to the end of the file
    sequences.sort(key=lambda record: node_to_date.get(record.id, float('inf')))

    with open(output_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python Sort_dates_3.py reference_file input_fasta input_dates.csv sorted_output.fasta")
        sys.exit(1)

    reference_file = sys.argv[1]
    fasta_file = sys.argv[2]
    csv_file = sys.argv[3]
    output_file = sys.argv[4]

    temp_sequences = add_sequences(reference_file, fasta_file)
    sort_sequences_by_date(temp_sequences, csv_file, output_file)
    print("Sorting Completed")
