import os
import glob
import csv
from collections import defaultdict

# Directory with besthit BED files
input_dir = 'stage_04_1'

# Output TSV file for besthit scores
output_file = 'gene_besthit_matrix_max.tsv'

# Nested dictionary: {gene_id: {sample_name: max_score}}
gene_matrix = defaultdict(dict)
all_gene_ids = set()
all_samples = []

# Find all besthit files (adjust pattern as needed)
input_files = glob.glob(os.path.join(input_dir, 'besthit-logpval@*@promoters_250u_10d_1.bed'))

for filepath in input_files:
    filename = os.path.basename(filepath)

    # Extract sample name from filename
    sample_name = filename.split('@')[1]
    all_samples.append(sample_name)

    gene_to_scores = defaultdict(list)

    with open(filepath, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue
            gene_id = cols[3]
            try:
                score = float(cols[4])
            except ValueError:
                continue
            gene_to_scores[gene_id].append(score)

    # Take maximum score per gene for this sample
    for gene_id, scores in gene_to_scores.items():
        max_score = max(scores)
        gene_matrix[gene_id][sample_name] = max_score
        all_gene_ids.add(gene_id)

# Write the gene besthit matrix to TSV
with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['gene_id'] + all_samples)

    for gene_id in all_gene_ids:
        row = [gene_id]
        for sample in all_samples:
            value = gene_matrix[gene_id].get(sample, 0.0)
            row.append(f"{value:.6E}")
        writer.writerow(row)

print(f"Done: saved gene besthit matrix to {output_file}")
