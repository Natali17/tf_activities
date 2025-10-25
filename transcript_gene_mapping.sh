#!/bin/bash
set -euo pipefail

GTF_FILE="../gencode.v44.annotation.gtf"
OUTPUT_FILE="transcript_to_gene.tsv"

echo -e "transcript_id\tgene_id" > "$OUTPUT_FILE"

# Deriving transcript_id и gene_id
awk 'BEGIN{OFS="\t"} 
     !/^#/ {
         if($9 ~ /transcript_id/ && $9 ~ /gene_id/) {
             match($9, /transcript_id "([^"]+)"/, tid)
             match($9, /gene_id "([^"]+)"/, gid)
             if(tid[1]!="" && gid[1]!="") {
                 split(tid[1], tparts, "\\.")
                 split(gid[1], gparts, "\\.")
                 print tparts[1], gparts[1]
             }
         }
     }' "$GTF_FILE" >> "$OUTPUT_FILE"

echo "Saved $(($(wc -l < "$OUTPUT_FILE") - 1)) correspondences in $OUTPUT_FILE"
