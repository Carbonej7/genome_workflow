import sys
import pandas as pd

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

# Retrieve input and output file names from command-line arguments
input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

# Read data from CSV with column headers
data = pd.read_csv(input_file_name)

# Extract only the first 4 digits from the numeric date column as 'Year'
data['Year'] = data['numeric date'].astype(str).str[:4]

# Get the column names for calculation (all columns starting from the 2nd column)
value_column_names = data.columns[1:]

# Calculate the average values for each year
unique_years = data['Year'].unique()

# Initialize a dictionary to store the individual average values
individual_average_values = {}

for year in unique_years:
    idx = data['Year'] == year
    individual_average_values[year] = data.loc[idx, value_column_names].mean()

# Create a DataFrame for the output
output_df = pd.DataFrame(individual_average_values).T.reset_index()

# Rename the index column to 'Year'
output_df = output_df.rename(columns={'index': 'Year'})

# Write the DataFrame to a CSV file
output_df.to_csv(output_file_name, index=False)

print(f'Individual average values per Year saved to "{output_file_name}"')
