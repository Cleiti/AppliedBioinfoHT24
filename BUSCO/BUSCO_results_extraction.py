import os
import re
import pandas as pd


filtered_dir = "/proj/uppstore2019013/nobackup/private/1000spider_master_project/busco/busco_results_filtered/res"
unfiltered_dir = "/proj/uppstore2019013/nobackup/private/1000spider_master_project/busco/busco_results_unfiltered"
output_csv = "/proj/uppstore2019013/nobackup/private/1000spider_master_project/busco/busco_comparison_results.csv"

#function to parse species name from filename
def extract_species_name(filename):
    #Remove prefixes like 'filtered_' and suffixes like '_out'
    return re.sub(r"^filtered_|_out$", "", filename)

# Helper function to extract a percentage from a line
def extract_percentage(pattern, text):
    match = re.search(pattern, text)
    return float(match.group(1)) if match else None

# Initialize a list to store data
busco_data = []

# Iterate over files in the filtered directory
for filtered_file in os.listdir(filtered_dir):
    if filtered_file.endswith("_out"):
        species = extract_species_name(filtered_file)

        # Determine the corresponding unfiltered file
        unfiltered_file = os.path.join(unfiltered_dir, f"{species}_out")
        filtered_file_path = os.path.join(filtered_dir, filtered_file)

        # Check if both filtered and unfiltered files exist
        if os.path.exists(unfiltered_file):
            with open(filtered_file_path, "r") as f:
                filtered_lines = f.read()
            with open(unfiltered_file, "r") as f:
                unfiltered_lines = f.read()

            # Extract percentages and values from both files
            filtered_completeness = extract_percentage(r"\|C:([\d.]+)%", filtered_lines)
            filtered_single_copy = extract_percentage(r"\[S:([\d.]+)%", filtered_lines)
            unfiltered_completeness = extract_percentage(r"\|C:([\d.]+)%", unfiltered_lines)
            unfiltered_single_copy = extract_percentage(r"\[S:([\d.]+)%", unfiltered_lines)

            # Handle missing values
            if None in {filtered_completeness, filtered_single_copy, unfiltered_completeness, unfiltered_single_copy}:
                print(f"Warning: Missing data for species '{species}'. Skipping...")
                continue

            # Add data to the list
            busco_data.append({
                "Species": species,
                "Completeness_unfiltered": unfiltered_completeness,
                "Completeness_filtered": filtered_completeness,
                "Single_copy_unfiltered": unfiltered_single_copy,
                "Single_copy_filtered": filtered_single_copy,
            })


df = pd.DataFrame(busco_data)
df.to_csv(output_csv, index=False)

print(f"BUSCO results comparison saved to {output_csv}")