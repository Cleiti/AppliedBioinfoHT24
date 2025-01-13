def filter_blast_results(input_file, output_file):
    seen = set() 
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            query_id = line.split()[0] 
            if query_id not in seen:
                outfile.write(line) 
                seen.add(query_id)

input_file = "reciprocal_clavata_4676.txt" 
output_file = "filtered_blast_results.txt" 

filter_blast_results(input_file, output_file)
