import re
import pandas as pd
import json

bigg_ref = pd.read_csv('bigg_models_reactions.txt', sep = '\t')
strains = ['Sq4a', 'GAOx1']

for strain in strains:
    blast_draft_in = f'../BLAST/draft_models/{strain}_draft_rxn_list.tsv'
    orphan_eggnog_in = f'eggnog_annotations/{strain}_eggnog_orphans_only.tsv'
    orphan_gene_set_in = f'orphan_reactions/{strain}_gene_set.txt'
    orphan_reaction_set_in = f'orphan_reactions/{strain}_reaction_set.txt'
    orphan_genes_reactions_in = f'orphan_reactions/{strain}_genes_reactions.json'

    orphans_only_draft_out = f'orphan_draft_models/{strain}_orphan_draft_rxn_list.tsv'
    combined_draft_out = f'../draft_reconstructions/draft_models_w_gaps/{strain}_draft_model_w_gaps.tsv'

    blast_draft = pd.read_csv(blast_draft_in, sep='\t')
    orphan_eggnog = pd.read_csv(orphan_eggnog_in, sep='\t')
    orphan_gene_set = []
    with open(orphan_gene_set_in, 'r') as file:
        for line in file:
            orphan_gene_set.append(line.rstrip())
    orphan_reaction_set = []
    with open(orphan_reaction_set_in, 'r') as file:
        for line in file:
            orphan_reaction_set.append(line.rstrip())
    with open(orphan_genes_reactions_in, 'r') as file:
        orphan_genes_reactions = json.load(file)
    
    blast_draft.columns = ['rxn_id', 'name', 'formula', 'gpr', 'gr', 'pr', 'subsystem', 'EC', 'rev', 'lower', 'upper', 'obj']
    new_reactions_only = []
    for rxn in orphan_reaction_set:
        if rxn in blast_draft['rxn_id']:
            continue
        new_reactions_only.append(rxn)
    new_draft_rxn_list = pd.DataFrame(columns=['rxn_id', 'name', 'formula', 'gr', 'BRITE_terms', 'EC', 'rev', 'lower', 'upper'])
    for rxn in new_reactions_only:
        bigg_ref_rxn = bigg_ref[bigg_ref['bigg_id'] == rxn]
        if len(bigg_ref_rxn) > 1:
            print("ERROR MORE THAN ONE MATCH")
            continue
        elif len(bigg_ref_rxn) == 0:
            bigg_ref_rxn = bigg_ref[bigg_ref['old_bigg_ids'].str.contains(rxn, na=False)]
            if len(bigg_ref_rxn) == 0:
                print("ERROR NO REACTION FOUND: ", rxn)
                continue
            rxn = bigg_ref_rxn['bigg_id']
        rxn_id = rxn
        rxn_name = bigg_ref_rxn['name'].iloc[0]
        rxn_formula = bigg_ref_rxn['reaction_string'].iloc[0]
        rxn_gr = orphan_genes_reactions[rxn]
        rxn_brite_terms = set()
        rxn_ec_num = set()
        for gene in rxn_gr:
            brite_terms = orphan_eggnog.loc[orphan_eggnog['locus_tag'] == gene, 'BRITE'].iloc[0].split(',')
            rxn_brite_terms.update(brite_terms)
            ec_num = orphan_eggnog.loc[orphan_eggnog['locus_tag'] == gene, 'EC'].iloc[0]
            if ec_num != '-':
                rxn_ec_num.update(ec_num)
        rxn_brite_terms = list(rxn_brite_terms)
        rxn_ec_num = list(rxn_ec_num)
        rev = '-'
        upper = '-'
        lower = '-'
        new_row = [rxn_id, rxn_name, rxn_formula, rxn_gr, rxn_brite_terms, rxn_ec_num, rev, lower, upper]
        new_draft_rxn_list.loc[len(new_draft_rxn_list)] = new_row

    new_draft_rxn_list.to_csv(orphans_only_draft_out, sep='\t', index=False)

    blast_draft_new = blast_draft[['rxn_id', 'name', 'formula', 'gr', 'subsystem', 'EC', 'rev', 'lower', 'upper']]
    new_draft_rxn_list.columns = ['rxn_id', 'name', 'formula', 'gr', 'subsystem', 'EC', 'rev', 'lower', 'upper']

    combined_draft = pd.concat([blast_draft_new, new_draft_rxn_list])
    combined_draft.to_csv(combined_draft_out, sep='\t', index=False)

    gtf_gaps_only_in = f'../BLAST/draft_models/gap_genes/{strain}_gap_genes.tsv'
    eggnog_gaps_only_in = f'eggnog_annotations/{strain}_eggnog_orphans_only.tsv'

    all_gaps_out = f'../gaps/{strain}_gap_genes.tsv'
    transport_only_out = f'../gaps/transporter_genes/{strain}_transporters.tsv'
    non_transport_out = f'../gaps/non_transporters/{strain}_non_transporters.tsv'

    gtf_gaps_only = pd.read_csv(gtf_gaps_only_in, sep='\t')
    remaining_gaps = gtf_gaps_only[['locus_tag', 'product', 'protein_id', 'go_function', 'go_process', 'go_component', 'gene']]

    eggnog_gaps_only = pd.read_csv(eggnog_gaps_only_in, sep='\t')
    eggnog_gaps_only = eggnog_gaps_only[['locus_tag','GOs', 'Description', 'KEGG_ko', 'KEGG_Module', 'BRITE', 'KEGG_TC', 'CAZy', 'PFAMs']]

    gap_filling_df = pd.merge(remaining_gaps, eggnog_gaps_only, on='locus_tag')
    gap_filling_df = gap_filling_df[~gap_filling_df['locus_tag'].isin(orphan_gene_set)]

    transport_genes = gap_filling_df.query('KEGG_TC != "-"')
    non_transport_genes = gap_filling_df.query('KEGG_TC == "-"')

    gap_filling_df.to_csv(all_gaps_out, sep='\t', index=False)
    transport_genes.to_csv(transport_only_out, sep='\t', index=False)
    non_transport_genes.to_csv(non_transport_out, sep='\t', index=False)