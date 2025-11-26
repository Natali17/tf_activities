import argparse
import gzip
import os

def parse_gtf_attributes(attr_string):
    """
    Parses the 9th column (attributes string) of a GTF file into a dictionary.
    This helper function is robust for extracting attributes like transcript_name.
    """
    attrs = {}
    # Splits the string by semicolon (;) and strips quotes and whitespace
    for attr in attr_string.strip(';').split(';'):
        attr = attr.strip()
        if not attr:
            continue
        try:
            # Splits the key-value pair, ensuring only the first space separates them
            key, value = attr.split(None, 1)
            # Removes surrounding quotes from the value
            attrs[key] = value.strip('"')
        except ValueError:
            # Skip malformed attributes
            continue 
    return attrs

def create_id_mapping(gtf_file, output_file):
    """
    Creates the default mapping: Transcript ID (ENST...) -> Gene ID (ENSG...).
    Removes version suffixes (e.g., .1, .2) from IDs.
    """
    print(f"Creating Transcript ID -> Gene ID mapping. Output file: {output_file}")
    transcript_to_gene = {}
    # Handle gzipped GTF files
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    
    # Open and process the GTF file
    with open_func(gtf_file, 'rt') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            attr = fields[8]
            
            # Check for both ID attributes in the line
            if 'transcript_id' in attr and 'gene_id' in attr:
                try:
                    # Original logic using string splitting to extract and remove versions
                    transcript_id = attr.split('transcript_id "')[1].split('"')[0].split('.')[0]
                    gene_id = attr.split('gene_id "')[1].split('"')[0].split('.')[0]
                    
                    transcript_to_gene[transcript_id] = gene_id
                except IndexError:
                    # Skip lines with incomplete attributes
                    continue
    
    # Save the mapping to a TSV file
    with open(output_file, 'w') as out:
        out.write("transcript_id\tgene_id\n")
        for transcript_id, gene_id in transcript_to_gene.items():
            out.write(f"{transcript_id}\t{gene_id}\n")

    print(f"Successfully saved {len(transcript_to_gene)} ID mappings.")


def create_name_mapping(gtf_file, output_file):
    """
    Creates the Transcript Name (DDX11L2-202) -> Gene Name (DDX11L2) mapping.
    This requires using the 'transcript' feature and parsing attributes carefully.
    """
    print(f"Creating Transcript Name -> Gene Name mapping. Output file: {output_file}")
    transcript_name_to_gene_name = {}
    open_func = gzip.open if gtf_file.endswith('.gz') else open

    with open_func(gtf_file, 'rt') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature = fields[2]
            
            # Only process 'transcript' features, as they contain 'transcript_name'
            if feature == 'transcript':
                attrs = parse_gtf_attributes(fields[8])
                
                # Check for the required 'name' attributes
                if 'transcript_name' in attrs and 'gene_name' in attrs: 
                    transcript_name = attrs['transcript_name']
                    gene_name = attrs['gene_name']
                    
                    # Store the mapping (no version removal needed for names)
                    transcript_name_to_gene_name[transcript_name] = gene_name
    
    # Save the mapping to a TSV file
    with open(output_file, 'w') as out:
        out.write("transcript_name\tgene_name\n")
        for tr_name, g_name in transcript_name_to_gene_name.items():
            out.write(f"{tr_name}\t{g_name}\n")

    print(f"Successfully saved {len(transcript_name_to_gene_name)} Name mappings.")


def main():
    """
    Main function to configure argument parsing and select the mapping mode.
    """
    parser = argparse.ArgumentParser(
        description="Creates a transcript-to-gene mapping file from a GTF annotation file."
    )
    parser.add_argument("gtf_file", help="Input GTF file path (.gtf or .gtf.gz)")
    parser.add_argument("-o", "--output", default="transcript_to_gene.tsv", 
                        help="Output TSV filename. Default is 'transcript_to_gene.tsv'")
    parser.add_argument("--use-names", action="store_true",
                        help="If set, creates Transcript Name -> Gene Name mapping (e.g., DDX11L2-202 -> DDX11L2). Otherwise, creates ID -> ID mapping (ENST... -> ENSG...).")
    
    args = parser.parse_args()

    if args.use_names:
        # If the user specified --use-names, run the name mapping function
        if args.output == "transcript_to_gene.tsv":
            # Change default output name for clarity if the user didn't specify one
            args.output = "transcript_name_to_gene_name.tsv"
        create_name_mapping(args.gtf_file, args.output)
    else:
        # Default: run the ID mapping function
        create_id_mapping(args.gtf_file, args.output)

if __name__ == "__main__":
    main()

