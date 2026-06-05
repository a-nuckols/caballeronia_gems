import pandas as pd
import re
import json

def extract_protein_id(query):
    if not isinstance(query, str):
        return False
    match = re.search(r'WP_\d{9}\.\d+', query)
    if match:
        return match.group()
    else:
        return None

def extract_locus_tag(query, strain):
    if not isinstance(query, str):
        return False
    raw = ''
    if strain == 'Sq4a':
        raw = r'LDZ24_RS\d+'
    elif strain == 'GAOx1':
        raw = r'NCS66_RS\d+'
    match = re.search(raw, query)
    if match:
        return match.group()
    else:
        return None

strains = ['Sq4a', 'GAOx1']
xref = 'reac_xref.tsv'
xref_df = pd.read_csv(xref, sep='\t', comment='#', header=None, names=['source', 'ID', 'description'])
xref_df = xref_df.iloc[1:]

for strain in strains:
    print(f'{strain}:')
    eggnog_csv = f'eggnog_annotations/{strain}_eggnog.csv'
    metadata_tsv = f'../BLAST/genome_fasta/{strain}_metadata.tsv'
    orphans_tsv = f'../BLAST/draft_models/gap_genes/{strain}_gap_genes.tsv'

    eggnog = pd.read_csv(eggnog_csv)
    metadata = pd.read_csv(metadata_tsv, sep='\t')
    orphans = pd.read_csv(orphans_tsv, sep='\t')
    eggnog['protein_id'] = eggnog['#query'].apply(extract_protein_id)
    eggnog_protein = eggnog[~eggnog['protein_id'].isna()]
    eggnog_locus = eggnog[eggnog['protein_id'].isna()].drop(columns='protein_id')
    eggnog_locus['locus_tag'] = eggnog_locus['#query'].apply(extract_protein_id)
    eggnog_locus = eggnog_locus[~eggnog_locus['locus_tag'].isna()]
    eggnog_protein = pd.merge(eggnog_protein, metadata, how='left', on='protein_id')
    eggnog_locus = pd.merge(eggnog_locus, metadata, how='left', on='locus_tag')
    eggnog_locus = eggnog_locus[eggnog_protein.columns]
    eggnog = pd.concat([eggnog_protein, eggnog_locus])
    eggnog = eggnog[eggnog['locus_tag'].isin(orphans['gene_id'])]
    
    eggnog.to_csv(f'eggnog_annotations/{strain}_eggnog_orphans_only.tsv', sep='\t', index=False)

    kegg_rxns = {}
    all_rxns = []

    for i in range(len(eggnog)):
        rxn_string = eggnog['KEGG_Reaction'].iloc[i]
        id = eggnog['locus_tag'].iloc[i]
        if rxn_string == '-':
            continue
        rxn_list = rxn_string.split(',')
        kegg_rxns[id] = rxn_list
        all_rxns.extend(set(rxn_list)-set(all_rxns))
    xref_df['source_db'] = ''
    db_pos = xref_df.columns.get_loc('source_db')
    xref_df['source_id'] = ''
    id_pos = xref_df.columns.get_loc('source_id')
    
    for i in range(len(xref_df)):
        source_str = xref_df.iloc[i, 0]
        source_split = source_str.split(':', 1)
        xref_df.iloc[i, db_pos] = source_split[0]
        if len(source_split) == 2:
            xref_df.iloc[i, id_pos] = source_split[1]
    
    kegg_xref = xref_df[xref_df["source_db"] == "kegg.reaction"]
    bigg_xref = xref_df[xref_df["source_db"] == "bigg.reaction"]
    biggR_xref = xref_df[xref_df["source_db"] == "biggR"]

    all_rxns_df = pd.DataFrame({'kegg_id': all_rxns})
    all_rxns_df = pd.merge(all_rxns_df, kegg_xref[['source_id', 'ID']], left_on = 'kegg_id', right_on = 'source_id').drop(columns='source_id')
    
    kegg_to_bigg = {}

    for i in range(len(all_rxns_df)):
        kegg_rxn = all_rxns_df['kegg_id'].iloc[i]
        mnx_id = all_rxns_df['ID'].iloc[i]
        bigg_rxns = bigg_xref.query('ID == @mnx_id')
        bigg_rxns = bigg_rxns[~(bigg_rxns['description']=='secondary/obsolete/fantasy identifier')]
        mask = bigg_rxns['description'].str.contains(
            r'biggC:(?:cx|cm|m|h|x|n|l|f|r|v|s|g)',
            na=False
        )
        bigg_rxns = bigg_rxns[~mask]
        if len(bigg_rxns) >= 1:
            kegg_to_bigg[kegg_rxn] = bigg_rxns.iloc[:, id_pos].to_list()
        elif len(bigg_rxns) == 0:
            biggR_rxns = biggR_xref.query('ID == @mnx_id')
            biggR_rxns = biggR_rxns[~(biggR_rxns['description']=='secondary/obsolete/fantasy identifier')]
            mask = biggR_rxns['description'].str.contains(
                r'biggC:(?:cx|cm|m|h|x|n|l|f|r|v|s|g)',
                na=False
            )
            biggR_rxns = biggR_rxns[~mask]
            if len(biggR_rxns) >= 1:
                kegg_to_bigg[kegg_rxn] = biggR_rxns.iloc[:, id_pos].to_list()
    bigg_rxn_list = []
    for rxn in all_rxns:
        try: 
            rxn_bigg = kegg_to_bigg[rxn]
        except KeyError:
            continue
        bigg_rxn_list.extend(set(rxn_bigg)-set(bigg_rxn_list))
    bigg_rxn_genes = {}
    for gene, kegg_list in kegg_rxns.items():
        for kegg_rxn in kegg_list:
            for bigg_rxn in kegg_to_bigg.get(kegg_rxn, []):
                bigg_rxn_genes.setdefault(bigg_rxn, []).append(gene)
    for rxn in bigg_rxn_genes:
        bigg_rxn_genes[rxn] = list(set(bigg_rxn_genes[rxn]))
    gene_set = set()
    for genes in bigg_rxn_genes.values():
        gene_set.update(genes)
    with open(f'orphan_reactions/{strain}_gene_set.txt', 'w') as file:
        file.write('\n'.join(gene_set))
    with open(f'orphan_reactions/{strain}_reaction_set.txt', 'w') as file:
        file.write('\n'.join(bigg_rxn_list))
    with open(f'orphan_reactions/{strain}_genes_reactions.json', 'w') as file:
        json.dump(bigg_rxn_genes, file, indent=4)
    print(f'NUMBER OF GENES: {len(gene_set)}')
    print(f'NUMBER OF REACTIONS: {len(bigg_rxn_list)}')          

                

        

    

    

