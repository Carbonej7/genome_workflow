import sys
import pandas as pd


def clean_dates(input_file, cleaned_output_file):
    # Read data into a DataFrame
    df = pd.read_csv(input_file, sep='\t')

    root_index = df['#node'].str.startswith('NODE-').idxmax()
    df.loc[root_index, '#node'] = 'root'

    df = df[~df['#node'].str.startswith('NODE_')]

    df.to_csv(cleaned_output_file, sep='\t', index=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_dates.tsv cleaned_output.tsv")
    else:
        input_file = sys.argv[1]
        cleaned_output_file = sys.argv[2]
        clean_dates(input_file, cleaned_output_file)
