import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Define function
def replace_insertions(ancestral_sequences_file, root_sequence_file, output_file):
    # Parse SeqIO from BioPython to read the defined files
    ancestral_sequences = SeqIO.to_dict(SeqIO.parse(ancestral_sequences_file, "fasta"))
    # Parse SeqIO that these are fasta files
    root_sequence = SeqIO.read(root_sequence_file, "fasta")

    # Define strings with SeqIO
    for key, ancestral_record in ancestral_sequences.items():
        ancestral_sequence = str(ancestral_record.seq)
        root_sequence_str = str(root_sequence.seq)

        if len(ancestral_sequence) != len(root_sequence_str):
            raise ValueError("Ancestral sequence and root sequence must have the same length")

        # Constructs new replaced_sequence with insertions from root_sequence:
        # zip() combines characters pairwise from str
        # root base if ancestral base == "-" or "N" else ancestral_base iterates over each base from pairwise analysis
        # and replaces them if they are "-" or "N" with the root base

        # .join() joins resulting sequences
        # Seq() creates Seq object from joined strings of characters converting to Biological sequence
        replaced_sequence = Seq("".join(root_base if ancestral_base == "-" or "N" else ancestral_base
                                        for ancestral_base, root_base in zip(ancestral_sequence, root_sequence_str)))

        ancestral_record.seq = replaced_sequence

    # Write new insertions file to output
    with open(output_file, "w") as output_handle:
        # Write root_sequence to output
        root_seq_record = SeqRecord(root_sequence.seq, id="root_sequence", description="root")
        SeqIO.write([root_seq_record], output_handle, "fasta")

        # Write modified ancestral sequences to output
        SeqIO.write(ancestral_sequences.values(), output_handle, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python replace_insertions.py ancestral_sequences.fasta root_sequence.fasta "
              "ancestral_seq_replaced_insertions.fasta")
        sys.exit(1)

    # Define inputs and outputs from command-line arguments
    ancestral_sequences_file = sys.argv[1]
    root_sequence_file = sys.argv[2]
    output_file = sys.argv[3]

    replace_insertions(ancestral_sequences_file, root_sequence_file, output_file)
