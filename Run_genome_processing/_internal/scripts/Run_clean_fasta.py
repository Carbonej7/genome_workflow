import re
import sys


def clean_accession_codes(fasta_file):
    cleaned_sequences = []

    with open(fasta_file, 'r') as file:
        lines = file.readlines()

        for line in lines:
            if line.startswith('>'):
                # Extract accession code without version number using regex
                accession_code = re.match(r'^>([^\s.]+)', line).group(1)
                cleaned_sequences.append(f'>{accession_code}\n')
            else:
                cleaned_sequences.append(line)

    output_file = 'sequences_aln_clean.fasta'
    with open(output_file, 'w') as output:
        output.writelines(cleaned_sequences)

    return output_file


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python script.py input_fasta_file")
        sys.exit(1)

    input_fasta_file = sys.argv[1]
    cleaned_fasta_file = clean_accession_codes(input_fasta_file)
    print(f'Cleaned version numbers from fasta save')
