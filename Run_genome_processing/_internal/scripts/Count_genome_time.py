import sys
from Bio import SeqIO
import pandas as pd


def count_differences(seq1, seq2):
    """Count the number of differences between two sequences."""
    return sum(1 for a, b in zip(seq1, seq2) if a != b)


def compare_with_root(root_sequence, sequences):
    """Compare each sequence with the root sequence and count differences."""
    differences = []
    for name, sequence in sequences.items():
        diff_count = count_differences(root_sequence, sequence)
        differences.append({'#node': name, 'differences': diff_count})
    return differences


def fasta_to_dataframe(fasta_file):
    """Read FASTA file and return a DataFrame with sequence names and sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    return pd.DataFrame.from_dict(sequences, orient='index', columns=['Sequence'])


def main(fasta_file, dates_csv, output_csv):
    # Read FASTA file
    sequences_df = fasta_to_dataframe(fasta_file)

    # Get the root sequence (first sequence)
    root_sequence = sequences_df.iloc[0]['Sequence']

    # Compare each sequence with the root sequence
    differences = compare_with_root(root_sequence, sequences_df['Sequence'])

    # Create DataFrame from differences
    differences_df = pd.DataFrame(differences)

    # Read dates from CSV
    dates_df = pd.read_csv(dates_csv)

    # Extract first 4 digits from the numeric date column
    dates_df['numeric date'] = dates_df['numeric date'].astype(str).str[:4]

    # Merge differences dataframe with dates dataframe
    merged_df = pd.merge(differences_df, dates_df, on='#node')

    # Calculate total differences per year
    differences_per_year = merged_df.groupby('numeric date')['differences'].sum()

    # Count number of sequences per year
    sequences_per_year = merged_df.groupby('numeric date').size()

    # Calculate average differences per year
    average_differences_per_year = differences_per_year / sequences_per_year

    # Create DataFrame for output
    output_df = pd.DataFrame({'numeric date': average_differences_per_year.index,
                              'average differences': average_differences_per_year.values})

    # Export average differences per year to CSV
    output_df.to_csv(output_csv, index=False)

    print(f'Results saved to "{output_csv}"')


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input.fasta dates.csv output.csv")
        sys.exit(1)
    fasta_file = sys.argv[1]
    dates_csv = sys.argv[2]
    output_csv = sys.argv[3]
    main(fasta_file, dates_csv, output_csv)
