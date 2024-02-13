import sys


# Function to clean dates data
def clean_dates(input_file, cleaned_output_file, root_output_file):
    # Variables to store information about the root entry
    root_entry = None

    # Read dates data
    with open(input_file, 'r') as input_file:
        lines = input_file.readlines()

        # Find the root entry
        for line in lines:
            if line.startswith('NODE_'):
                if root_entry is None:
                    root_entry = line.rstrip().split('\t')
                else:
                    break

        # Write the root entry to root_output_file
        if root_entry:
            with open(root_output_file, 'w') as root_file:
                root_file.write('\t'.join(root_entry) + '\n')

        # Write the cleaned entries to cleaned_output_file
        with open(cleaned_output_file, 'w') as cleaned_file:
            for line in lines:
                if not line.startswith('NODE_'):
                    cleaned_file.write(line)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input_dates.tsv cleaned_output.tsv root_output.tsv")
    else:
        input_file = sys.argv[1]
        cleaned_output_file = sys.argv[2]
        root_output_file = sys.argv[3]
        clean_dates(input_file, cleaned_output_file, root_output_file)
