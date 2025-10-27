#!/bin/bash
set -ex

# Step 1: Run the annotation script to process the GTF file
python annotation.py gencode.v44.annotation.gtf \
    --annotation-source gencode

# Step 2: Process motif data
bash process_motifs.sh

# Step 3: Map transcripts to genes
python transcript_gene_mapping.py

# Step 4: Convert transcript-level occupancy to gene-level occupancy
python trans_occ_to_gene_occ.py

# Step 5: Calculate gene occupancy averages
python gene_occupancy_avg.py

# Step 5: Calculate gene occupancy averages
python filtered_gene_id.py tcga_log2tmm.tsv gene_occupancy_matrix_avg.tsv


