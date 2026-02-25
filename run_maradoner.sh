#!/bin/bash

# List of distance values (d)
distances=(2000)

# List of file prefixes (besthit and occupancy)
file_types=("gene_occupancy")

# Main loop through file types
for type in "${file_types[@]}"; do
    # Nested loop through distances
    for d in "${distances[@]}"; do
        
        # Define the input filename
        # Example: gene_besthit_500u_250d_filtered.tsv
        INPUT_FILE="${type}_500u_${d}d_filtered.tsv"
        
        # Define a unique Project ID for each run
        # Shortens "besthit" to "bh" and "occupancy" to "occ" for the project name
        if [[ "$type" == *"besthit"* ]]; then
            PROJECT_ID="sclc_bh_500u_${d}d"
        else
            PROJECT_ID="sclc_occ_500u_${d}d"
        fi

        # Check if the input file actually exists before running
        if [[ -f "$INPUT_FILE" ]]; then
            echo "===================================================="
            echo "STARTING PROCESS FOR: $INPUT_FILE"
            echo "PROJECT ID: $PROJECT_ID"
            echo "===================================================="

            # 1. Create project
            echo "[1/4] Creating maradoner project..."
            maradoner create "$PROJECT_ID" debatched_subtype_without_stage_41samples_filtered.tsv "$INPUT_FILE" \
                --sample-groups data/subtype_41samples.json \
                --loading-transform esf

            # 2. Fit model
            echo "[2/4] Fitting the model..."
            maradoner fit "$PROJECT_ID"

            # 3. Predict
            echo "[3/4] Running predictions..."
            maradoner predict "$PROJECT_ID"

            # 4. Export results
            echo "[4/4] Exporting results to folder: $PROJECT_ID"
            maradoner export "$PROJECT_ID" output/subtype_without_stage_41samples/debatched_subtype_without_stage_41samples

            echo -e "SUCCESS: Completed processing for $PROJECT_ID\n"
        else
            echo "SKIPPING: File $INPUT_FILE not found."
        fi
        
    done
done

echo "All tasks completed."
