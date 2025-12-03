import os
import glob
import csv
from collections import defaultdict
from statistics import mean

# Directory containing occupancy BED files
input_dir = 'stage_04_1'

# Output TSV file
output_file = 'gene_occupancy_avg.tsv'

# Nested dictionary to store average occupancy per gene per sample
# Structure: {gene_id: {sample_name: occupancy_avg}}
gene_matrix = defaultdict(dict)

# Sets to track all unique genes and samples
all_gene_ids = set()
all_samples = []

# Find all input BED files matching the pattern
input_files = glob.glob(os.path.join(input_dir, 'occupancy@*@promoters_250u_10d_1.bed'))

for filepath in input_files:
    filename = os.path.basename(filepath)

    # Extract sample name from filename between "occupancy@" and "@promoters"
    sample_name = filename.split('@')[1]
    all_samples.append(sample_name)

    # Temporary dictionary: {gene_id: list of occupancy values for this sample}
    gene_to_occupancy = defaultdict(list)

    # Read each line from the BED file
    with open(filepath, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue  # skip incomplete lines
            gene_id = cols[3]
            try:
                occupancy = float(cols[4])
            except ValueError:
                continue  # skip non-numeric occupancy values

            gene_to_occupancy[gene_id].append(occupancy)

    # Compute average occupancy per gene and store in the global matrix
    for gene_id, occupancies in gene_to_occupancy.items():
        avg_occupancy = mean(occupancies)
        gene_matrix[gene_id][sample_name] = avg_occupancy
        all_gene_ids.add(gene_id)

# Write the aggregated gene occupancy matrix to a TSV file
with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['gene_id'] + all_samples)

    for gene_id in all_gene_ids:
        row = [gene_id]
        for sample in all_samples:
            value = gene_matrix[gene_id].get(sample, 0.0)
            row.append(f"{value:.6E}")  # scientific notation
        writer.writerow(row)

print(f"Done: saved gene occupancy matrix to {output_file}")
