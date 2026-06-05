import pandas as pd
from Bio import SeqIO

def parse_attributes(attr_string):
    attrs = {}
    for field in attr_string.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) == 2:
            key, value = parts
            attrs[key] = value.strip('"')
    return attrs

strains = ['Sq4a', 'GAOx1']
for strain in strains:
    gtf = pd.read_csv(
        f"../genomes/{strain}/annotation/{strain}_annotations.gtf",
        sep="\t",
        comment="#",
        header=None,
        names=[
            "seqid",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes"
        ]
    )

    attr_df = gtf["attributes"].apply(parse_attributes).apply(pd.Series)

    gtf = pd.concat(
        [gtf.drop(columns="attributes"), attr_df],
        axis=1
    )

    blast_hits = pd.read_csv(
        f"draft_models/{strain}_draft_blast_hits.tsv",
        sep='\t'
    )
    fasta_dict = {}
    for record in SeqIO.parse(f'genome_fasta/{strain}_proteins.faa', 'fasta'):
        fasta_dict[record.id] = record.seq

    gtf_gaps = gtf[~gtf['gene_id'].isin(set(blast_hits['qseqid']))]
    gtf_gaps = gtf_gaps[gtf_gaps['feature']=='CDS']
    gtf_gaps['sequence'] = gtf_gaps['gene_id'].map(fasta_dict)
    gtf_gaps.to_csv(
        f'draft_models/gap_genes/{strain}_gap_genes.tsv',
        sep='\t',
        index=False
    )
