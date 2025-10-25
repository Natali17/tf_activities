import sys

def load_gene_ids(filename, gene_col=0):
    """
    Load a set of gene IDs from a tab-delimited file.

    Parameters:
        filename (str): Path to the input file.
        gene_col (int): Column index containing gene IDs (default: 0).

    Returns:
        set: A set of gene IDs.
    """
    gene_ids = set()
    with open(filename, 'r') as f:
        header = next(f)
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) > gene_col:
                gene_ids.add(cols[gene_col])
    return gene_ids

def filter_file_by_genes(input_file, output_file, genes_to_keep, gene_col=0):
    """
    Filter a tab-delimited file to keep only rows with gene IDs in genes_to_keep.

    Parameters:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
        genes_to_keep (set): Set of gene IDs to retain.
        gene_col (int): Column index containing gene IDs (default: 0).
    """
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            cols = line.strip().split('\t')
            if len(cols) > gene_col and cols[gene_col] in genes_to_keep:
                fout.write(line)

# Check command-line arguments
if len(sys.argv) != 3:
    print("Usage: python filtered_gene_id.py <file1> <file2>")
    sys.exit(1)

# Assign input file paths from command-line arguments
file1 = sys.argv[1]
file2 = sys.argv[2]

# Load gene ID sets from both files
genes1 = load_gene_ids(file1)
genes2 = load_gene_ids(file2)

# Compute intersection of gene IDs
common_genes = genes1.intersection(genes2)
print(f"Number of common gene IDs: {len(common_genes)}")

# Filter both files based on common gene IDs and save output
filter_file_by_genes(file1, f"{file1.rsplit('.', 1)[0]}_filtered.tsv", common_genes)
filter_file_by_genes(file2, f"{file2.rsplit('.', 1)[0]}_filtered.tsv", common_genes)

print("Filtered files containing only common gene IDs have been saved.")
