import os
import glob
import argparse

def load_mapping(mapping_file, use_names):
    """
    Loads the appropriate mapping dictionary based on the --use-names flag.
    
    If use_names is True, it loads Transcript Name -> Gene Name.
    If use_names is False, it loads Transcript ID -> Gene ID (handling versions).
    """
    mapping = {}
    
    # Check if the mapping file exists
    if not os.path.exists(mapping_file):
        print(f"Error: Mapping file '{mapping_file}' not found.")
        print("Please run transcript_gene_mapping.py first to create it.")
        return mapping
        
    print(f"Loading mapping from {mapping_file}...")
    
    with open(mapping_file, 'r') as f:
        # Skip header line
        next(f) 
        
        for line in f:
            try:
                # Expecting two columns: [Transcript/Name] \t [Gene/Name]
                key, value = line.strip().split('\t')
                
                if not use_names:
                    # Original logic: ID -> ID mapping. Remove version suffix from the key.
                    base_key = key.split('.')[0] 
                    mapping[base_key] = value.split('.')[0] # Remove version from gene ID as well, though usually gene IDs are unversioned in the mapping file key
                else:
                    # Name -> Name mapping. Names do not contain version suffixes.
                    mapping[key] = value
            except ValueError:
                print(f"Skipping malformed line in mapping file: {line.strip()}")
                continue
                
    print(f"Loaded {len(mapping)} unique mappings.")
    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Converts transcript IDs/names in BED files to gene IDs/names."
    )
    parser.add_argument("--use-names", action="store_true",
                        help="If set, uses Transcript Name -> Gene Name mapping. Otherwise, uses Transcript ID -> Gene ID.")
    
    args = parser.parse_args()

    # Determine mapping file and key processing based on argument
    if args.use_names:
        mapping_file = 'transcript_name_to_gene_name.tsv'
        # In this mode, we expect the input ID to be the full Transcript Name (no version)
    else:
        mapping_file = 'transcript_to_gene.tsv'
        # In this mode, we expect the input ID to have a version suffix that needs removal

    transcript_to_gene_map = load_mapping(mapping_file, args.use_names)
    
    if not transcript_to_gene_map:
        return # Stop execution if mapping failed to load

    # Create output directory if it doesn't exist
    os.makedirs('stage_04_1', exist_ok=True)

    # Find all input files matching the pattern
    input_files = glob.glob('stages/stage_04/besthit-logpval@*@promoters_250u_10d.bed')

    if not input_files:
        print("Warning: No input files found in 'stage_04/'. Check the path and file pattern.")
        return

    for input_file in input_files:
        filename = os.path.basename(input_file)
        # Append '_1' to the filename for the output file
        output_filename = filename.replace('.bed', '_1.bed')
        output_path = os.path.join('stage_04_1', output_filename)

        with open(input_file, 'r') as fin, open(output_path, 'w') as fout:
            print(f"Processing {filename}...")
            lines_processed = 0
            lines_mapped = 0
            
            for line in fin:
                cols = line.strip().split('\t')
                lines_processed += 1
                
                if len(cols) >= 4:
                    input_id = cols[3]
                    
                    # --- Determine mapping key ---
                    if args.use_names:
                        # Key is the full Transcript Name (e.g., "DDX11L2-202")
                        map_key = input_id
                    else:
                        # Key is the Transcript ID base (e.g., "ENST00000456328")
                        map_key = input_id.split('.')[0] 
                    
                    # --- Perform mapping ---
                    # Replace input_id (col 3) with gene ID/Name; fallback to original if not found
                    gene_identifier = transcript_to_gene_map.get(map_key, input_id)
                    
                    if gene_identifier != input_id:
                        lines_mapped += 1
                        
                    cols[3] = gene_identifier
                
                fout.write('\t'.join(cols) + '\n')
            
            print(f"Finished. Total lines: {lines_processed}. Mapped lines: {lines_mapped}.")


if __name__ == "__main__":
    main()
