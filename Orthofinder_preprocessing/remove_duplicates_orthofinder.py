import pandas as pd
import re

# File paths
input_file = "/proj/uppstore2019013/nobackup/private/1000spider_master_project/orthofinder_filtered/filtered_with_ids.tsv"
output_file = "/proj/uppstore2019013/nobackup/private/1000spider_master_project/orthofinder_filtered/filtered_with_ids_no_dupes.tsv"

# Load the file
try:
    # Read the TSV file
    data = pd.read_csv(input_file, sep="\t", index_col=0)

    # Count initial number of columns
    initial_columns = data.shape[1]
    print(f"Initial number of columns: {initial_columns}")

    # Filter out columns ending with a period and a number
    filtered_data = data.loc[:, ~data.columns.str.contains(r"\.\d+$")]

    # Count remaining columns
    final_columns = filtered_data.shape[1]
    print(f"Number of columns after removing duplicates: {final_columns}")

    # Save the filtered data
    filtered_data.to_csv(output_file, sep="\t")
    print(f"Filtered file saved to {output_file}")
except FileNotFoundError:
    print(f"Error: File '{input_file}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")