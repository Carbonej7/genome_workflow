import sys
import pandas as pd


def tsv_to_csv(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    df.to_csv(output_file, index=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_tsv_to_csv.py <input_tsv_file> <output_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    tsv_to_csv(input_file, output_file)
