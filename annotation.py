import argparse
import gzip

class GTFRecord:
    """Represent a single GTF record."""
    def __init__(self, contig, source, feature, start, end, strand, attributes):
        self.contig = contig
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.attributes = attributes


def parse_gtf(gtf_file, tsl_filter=None, mane_filter=None, annotation_source='gencode'):
    """
    Parse a GTF file and return dictionaries of transcripts and genes.

    Parameters:
        gtf_file (str): Path to GTF (.gtf or .gtf.gz)
        tsl_filter (list of str): Optional transcript_support_level values to include
        mane_filter (bool): If True, include only transcripts with MANE Select tag
        annotation_source (str): Annotation source format

    Returns:
        transcripts (dict): transcript_id -> transcript info
        genes (dict): gene_id -> gene info
    """
    transcripts = {}
    genes = {}

    open_func = gzip.open if gtf_file.endswith('.gz') else open
    with open_func(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            attrs = {}
            for attr in fields[8].strip(';').split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                key, value = attr.split(None, 1)
                attrs[key] = value.strip('"')

            rec = GTFRecord(
                contig=fields[0],
                source=fields[1],
                feature=fields[2],
                start=fields[3],
                end=fields[4],
                strand=fields[6],
                attributes=attrs
            )

            # Apply TSL filter if provided
            if tsl_filter is not None and rec.feature == 'transcript':
                tsl = attrs.get('transcript_support_level')
                if tsl not in tsl_filter:
                    continue

            if mane_filter and 'MANE Select' not in attrs.get('tag', ''):
                continue

            if rec.feature == 'transcript':
                transcript_id = attrs.get('transcript_id')
                gene_id = attrs.get('gene_id')
                gene_name = attrs.get('gene_name', gene_id)
                if transcript_id and gene_id:
                    transcripts[transcript_id] = {
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'contig': rec.contig,
                        'start': rec.start,
                        'end': rec.end,
                        'strand': rec.strand,
                        'tsl': attrs.get('transcript_support_level', 'NA')
                    }

            elif rec.feature == 'gene':
                gene_id = attrs.get('gene_id')
                gene_name = attrs.get('gene_name', gene_id)
                if gene_id:
                    genes[gene_id] = {
                        'gene_name': gene_name,
                        'contig': rec.contig,
                        'start': rec.start,
                        'end': rec.end,
                        'strand': rec.strand
                    }

    return transcripts, genes


def get_tss(transcript):
    """Return transcription start site (TSS) based on strand."""
    return transcript['start'] if transcript['strand'] == '+' else transcript['end']


def configure_argparser():
    """Setup command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Extract TSS from GTF with optional TSL filtering and save as BED"
    )
    parser.add_argument("gtf_annotation", help="Input GTF file (.gtf or .gtf.gz)")
    parser.add_argument("-o", "--output", default="tss.bed", help="Output BED file")
    parser.add_argument("--tsl", type=int, nargs="+", default=None,
                        help="Filter by transcript_support_level (e.g., --tsl 1 2)")
    parser.add_argument("--filter-mane-select", action="store_true",
                        help="Filter transcripts to include only those with 'mane_select' tag (used "
                        "with MANE annotation)")
    parser.add_argument("--annotation-source", choices=["ensembl", "gencode", "mane"], default="gencode",
                        help="Specify the annotation source format for parsing (default: gencode)")
    parser.add_argument("--compress", action="store_true", help="Compress output as gzip (.gz)")
    return parser


def main():
    args = configure_argparser().parse_args()

    # Convert TSL filter to strings (GTF uses string values like "1", "2", etc.)
    tsl_filter = None
    if args.tsl is not None:
        tsl_filter = [str(tsl) for tsl in args.tsl]

    transcripts, genes = parse_gtf(args.gtf_annotation, tsl_filter=tsl_filter, annotation_source=args.annotation_source,
                            mane_filter=args.filter_mane_select)

    # Always output BED (6 columns)
    open_out = gzip.open if args.compress else open
    mode = 'wt' if args.compress else 'w'
    with open_out(args.output, mode) as f_out:
        for transcript_id, transcript in transcripts.items():
            tss = get_tss(transcript)
            contig = transcript['contig']
            strand = transcript['strand']
            chrom_start = tss - 1  # BED is 0-based, half-open
            chrom_end = tss
            name = transcript_id
            score = "0"
            if not contig.startswith("chr"):
                contig = f"chr{contig}"
            f_out.write(f"{contig}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\n")



if __name__ == "__main__":
    main()
