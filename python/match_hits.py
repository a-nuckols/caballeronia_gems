import pandas as pd
import numpy as np
import re

def eval_rule(rule, genes):
    if not isinstance(rule, str):
        return False
    genes = set(genes)
    tokens = set(re.findall(r'[a-z]\d+', rule))
    for gene in tokens:
        rule = re.sub(rf'\b{gene}\b', str(gene in genes), rule)
    try:
        return eval(rule)
    except Exception:
        return False

iJO1366 = pd.read_csv('iJO1366/iJO1366_rxn_list.tsv', sep="\t")

strains = ['Sq4a', 'GAOx1']

for strain in strains:
    as_query_tsv = f"blast_results/{strain}/filtered/{strain}_against_MG1655_filtered.tsv"
    as_target_tsv = f"blast_results/{strain}/filtered/MG1655_against_{strain}_filtered.tsv"
    draft_output_tsv = f"draft_models/{strain}_draft_rxn_list.tsv"
    final_blast_out = f"draft_models/{strain}_draft_blast_hits.tsv"
    final_blast_out2 = f"blast_results/{strain}/matched/{strain}_draft_blast_hits.tsv"

    as_query = pd.read_csv(as_query_tsv, sep='\t')
    as_target = pd.read_csv(as_target_tsv, sep='\t')

    final_blast_df = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    gene_list = []
    gene_dict = {}

    for id in pd.unique(as_query['sseqid']):
        as_query_hits = as_query.query(f'sseqid==@id')
        as_target_hits = as_target.query(f'qseqid==@id')
        if (len(as_query_hits)==0 or len(as_target_hits)==0):
            continue
        matches = as_query_hits[as_query_hits['qseqid'].isin(as_target_hits['sseqid'])]
        final_blast_df = pd.concat([final_blast_df, matches], ignore_index=True)
        if len(matches) > 0:
            gene_list.append(id)
            gene_dict[id] = matches['qseqid'].tolist()
    
    gene_set = set(gene_list)

    mask = iJO1366['Gene-Reaction Association'].apply(
        lambda rule: eval_rule(rule, gene_set) if isinstance(rule, str) else False
    )

    strain_draft_model = iJO1366[mask]
    print(f'{strain}: ')
    print(f'total reactions kept: {len(strain_draft_model)}')    
    strain_draft_model.to_csv(
        draft_output_tsv,
        sep='\t',
        index=False
    )
    final_blast_df.to_csv(
        final_blast_out,
        sep='\t',
        index=False
    )
    final_blast_df.to_csv(
        final_blast_out2,
        sep='\t',
        index=False
    )



