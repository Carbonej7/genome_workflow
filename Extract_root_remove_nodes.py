import sys
from Bio import SeqIO

# Function to export fasta data to a file
def export_fasta(filename, header, sequence):
    with open(filename, 'w') as file:
        file.write(f'>{header}\n{sequence}\n')

# Function to clean fasta data
def clean_fasta(input_file, cleaned_output_file, root_output_file):
    # Variables to store information about the root sequence
    root_header = None
    root_sequence = ""

    # Read fasta data
    for seq_record in SeqIO.parse(input_file, "fasta"):
        if seq_record.id.startswith('NODE_') and root_header is None:
            root_header = 'root'
            root_sequence = str(seq_record.seq)
            break

    # Export root sequence to root.fasta
    if root_header is not None and root_sequence:
        export_fasta(root_output_file, root_header, root_sequence)

    # Remove root sequence from the rest of the fasta data and write cleaned data to a file
    with open(cleaned_output_file, 'w') as cleaned_file:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            if not seq_record.id.startswith('NODE_'):
                cleaned_file.write(f">{seq_record.id}\n{seq_record.seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input.fasta cleaned_output.fasta root_output.fasta")
    else:
        input_file = sys.argv[1]
        cleaned_output_file = sys.argv[2]
        root_output_file = sys.argv[3]
        clean_fasta(input_file, cleaned_output_file, root_output_file)

