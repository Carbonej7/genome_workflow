import sys

def read_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as f:
        header = None
        sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    sequences[header] = sequence
                header = line[1:]
                sequence = ''
            else:
                sequence += line
        if header is not None:
            sequences[header] = sequence
    return sequences

def read_tsv(tsv_file):
    sequences_to_remove = set()
    with open(tsv_file, 'r') as f:
        for line in f:
            sequence = line.strip()
            sequences_to_remove.add(sequence)
    return sequences_to_remove

def write_fasta(fasta_file, sequences):
    with open(fasta_file, 'w') as f:
        for header, sequence in sequences.items():
            f.write('>' + header + '\n')
            f.write(sequence + '\n')

def remove_sequences(fasta_file, tsv_file, output_file):
    fasta_sequences = read_fasta(fasta_file)
    sequences_to_remove = read_tsv(tsv_file)

    filtered_sequences = {header: sequence for header, sequence in fasta_sequences.items() if header not in sequences_to_remove}

    write_fasta(output_file, filtered_sequences)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input.fasta sequences_to_remove.tsv output.fasta")
    else:
        input_fasta = sys.argv[1]
        sequences_to_remove_tsv = sys.argv[2]
        output_fasta = sys.argv[3]
        remove_sequences(input_fasta, sequences_to_remove_tsv, output_fasta)
