gtf_file = "gencode.v44.annotation.gtf"
output_file = "transcript_to_gene.tsv"

transcript_to_gene = {}

with open(gtf_file, 'r') as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        attr = fields[8]
        if 'transcript_id' in attr and 'gene_id' in attr:
            try:
                transcript_id = attr.split('transcript_id "')[1].split('"')[0].split('.')[0]
                gene_id = attr.split('gene_id "')[1].split('"')[0].split('.')[0]
                transcript_to_gene[transcript_id] = gene_id
            except IndexError:
                continue

# Сохраняем в .tsv
with open(output_file, 'w') as out:
    out.write("transcript_id\tgene_id\n")
    for transcript_id, gene_id in transcript_to_gene.items():
        out.write(f"{transcript_id}\t{gene_id}\n")

print(f"Сохранено {len(transcript_to_gene)} соответствий в файл {output_file}")

