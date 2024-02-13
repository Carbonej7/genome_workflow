import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def replace_insertions(root_file, ref_sequence_file, output_file):
    root_records = SeqIO.to_dict(SeqIO.parse(root_file, "fasta"))
    ref_sequence = SeqIO.read(ref_sequence_file, "fasta")

    for key, root_record in root_records.items():
        root_seq = str(root_record.seq)
        ref_sequence_str = str(ref_sequence.seq)

        # Replace insertions in ancestral sequence with corresponding points from ref sequence
        replaced_sequence = Seq("".join(ref_base if root_base == "-" or "N" else root_base
                                        for root_base, ref_base in zip(root_seq, ref_sequence_str)))

        root_record.seq = replaced_sequence

    # Write the modified ancestral sequences to the output file
    SeqIO.write(list(root_records.values()), output_file, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python replace_insertions.py ancestral_sequences.fasta ref_sequence.fasta "
              "ancestral_seq_replaced_insertions.fasta")
        sys.exit(1)

    root_file = sys.argv[1]
    ref_sequence_file = sys.argv[2]
    output_file = sys.argv[3]

    replace_insertions(root_file, ref_sequence_file, output_file)
