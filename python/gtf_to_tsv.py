import pandas as pd

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

    gtf_gaps = gtf[~gtf['gene_id'].isin(set(blast_hits['qseqid']))]

    gtf_gaps.to_csv(
        f'draft_models/gap_genes/{strain}_gap_genes.tsv',
        sep='\t',
        index=False
    )
