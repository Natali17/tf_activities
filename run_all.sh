#!/bin/bash
set -ex

# Step 1: Run the annotation script to process the GTF file
python gtf_mapping.py gencode.v44.annotation.gtf \
    --annotation-source gencode --use-names

# Step 2: Process motif data
bash process_motifs.sh # {pipeline_core.sh, make_flanks.sh}

# Step 3: Map transcripts to genes
python transcript_gene_mapping.py 

# Step 4: Convert transcript-level occupancy to gene-level occupancy
python trans_occ_to_gene_occ.py --use-names

# Step 5: Calculate gene occupancy averages
python aggregate_by_metrics.py #OR aggregate_less_memory.py for big .tsv files

# Step 5: Calculate gene occupancy averages
bash filtering_pipeline.sh #


