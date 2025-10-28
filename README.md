# Description
This repository contains scripts and analysis pipelines for transcriptional motif analysis and gene expression normalization (TMM-based preprocessing) with the use of [MARADONER](motif activities)
Script run_all.sh helps to understand the sequence of steps to preprocess the data:

## Contents
### Annotation processing:
- `annotation.py`: annotation of GTF files
**You can download GENCODE annotation data** 
```python
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gzip -d gencode.v44.annotation.gtf.gz
```

- `process_motifs.sh`: processing motif data
- `transcript_gene_mapping.py`: mapping transcripts to genes
- `trans_occ_to_gene_occ.py`: transcript-to-gene occupancy transformation
- `gene_occupancy_avg.py`: gene occupancy averages calculation
  OR `gene_besthit_max.py`: gene besthit maximum calculation

### Expression data processing
- `tcga_counts_to_log2tmm.R`: TCGA expression preprocessing
- `cell_lines_counts_to_log2tmm.R`: cell lines expression preprocessing

### Data filtering to avoid errors during maradoner running:
- `filtered_gene_id.py` <file1> <file2> OR `filtered_gene_id_3files.py` <file1> <file2> <file3>: leaving only those gene_ids that match in the specified files

All the scripts should be in the same directory to work correctly
