import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def replace_insertions(ancestral_sequences_file, root_sequence_file, output_file):
    ancestral_sequences = SeqIO.to_dict(SeqIO.parse(ancestral_sequences_file, "fasta"))
    root_sequence = SeqIO.read(root_sequence_file, "fasta")

    for key, ancestral_record in ancestral_sequences.items():
        ancestral_sequence = str(ancestral_record.seq)
        root_sequence_str = str(root_sequence.seq)

        if len(ancestral_sequence) != len(root_sequence_str):
            raise ValueError("Ancestral sequence and root sequence must have the same length")

        # Replace insertions in ancestral sequence with corresponding bases from root sequence
        replaced_sequence = Seq("".join(root_base if ancestral_base == "-" else ancestral_base
                                        for ancestral_base, root_base in zip(ancestral_sequence, root_sequence_str)))

        ancestral_record.seq = replaced_sequence

    # Write the modified ancestral sequences to the output file
    SeqIO.write(list(ancestral_sequences.values()), output_file, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python replace_insertions.py ancestral_sequences.fasta root_sequence.fasta "
              "ancestral_seq_replaced_insertions.fasta")
        sys.exit(1)

    ancestral_sequences_file = sys.argv[1]
    root_sequence_file = sys.argv[2]
    output_file = sys.argv[3]

    replace_insertions(ancestral_sequences_file, root_sequence_file, output_file)
