import os
import glob

# Load transcript-to-gene mapping without versions
transcript_to_gene = {}
with open('transcript_to_gene.tsv', 'r') as f:
    for line in f:
        transcript, gene = line.strip().split('\t')
        base_transcript = transcript.split('.')[0]  # remove version suffix
        transcript_to_gene[base_transcript] = gene

# Create output directory if it doesn't exist
os.makedirs('stage_04_1', exist_ok=True)

# Find all input files matching the pattern
input_files = glob.glob('../stage_04/occupancy@*@promoters_250u_10d.bed')

for input_file in input_files:
    filename = os.path.basename(input_file)
    output_filename = filename.replace('.bed', '_1.bed')
    output_path = os.path.join('stage_04_1', output_filename)

    with open(input_file, 'r') as fin, open(output_path, 'w') as fout:
        for line in fin:
            cols = line.strip().split('\t')
            if len(cols) >= 4:
                full_transcript_id = cols[3]
                base_transcript_id = full_transcript_id.split('.')[0]  # remove version suffix
                # Replace transcript ID with corresponding gene ID; fallback to original transcript ID
                gene_id = transcript_to_gene.get(base_transcript_id, full_transcript_id)
                cols[3] = gene_id
            fout.write('\t'.join(cols) + '\n')
