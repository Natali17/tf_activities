# Define input GTF file and output mapping file
gtf_file = "gencode.v44.annotation.gtf"
output_file = "transcript_to_gene.tsv"

# Dictionary to store transcript-to-gene mappings
transcript_to_gene = {}

# Open and parse the GTF annotation file
with open(gtf_file, 'r') as gtf:
    for line in gtf:
        # Skip comment lines
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        # Skip malformed lines
        if len(fields) < 9:
            continue
        
        attr = fields[8]
        # Extract transcript_id and gene_id attributes
        if 'transcript_id' in attr and 'gene_id' in attr:
            try:
                # Remove version numbers (e.g., .1, .2)
                transcript_id = attr.split('transcript_id "')[1].split('"')[0].split('.')[0]
                gene_id = attr.split('gene_id "')[1].split('"')[0].split('.')[0]
                transcript_to_gene[transcript_id] = gene_id
            except IndexError:
                # Skip lines with incomplete attributes
                continue

# Save the transcript-to-gene mapping to a TSV file
with open(output_file, 'w') as out:
    out.write("transcript_id\tgene_id\n")
    for transcript_id, gene_id in transcript_to_gene.items():
        out.write(f"{transcript_id}\t{gene_id}\n")

# Print summary information
print(f"Saved {len(transcript_to_gene)} transcript-to-gene mappings to {output_file}")

