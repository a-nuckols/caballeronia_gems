import pandas as pd
import numpy as np
from Bio import SeqIO

strains = ['Sq4a', 'GAOx1']
MG1655_fasta = "genome_fasta/MG1655_proteins.faa"
MG1655_dict = {}
for record in SeqIO.parse(MG1655_fasta, 'fasta'):
    id = record.name
    value = len(record.seq)
    MG1655_dict[id]=value

for strain in strains:
    print(f"{strain}:")
    as_query_tsv = f'blast_results/{strain}/{strain}_against_MG1655'
    as_target_tsv = f'blast_results/{strain}/MG1655_against_{strain}'

    as_query_df = pd.read_csv(as_query_tsv, sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    as_target_df = pd.read_csv(as_target_tsv, sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    
    strain_fasta = f'genome_fasta/{strain}_proteins.faa'
    strain_dict = {}

    for record in SeqIO.parse(strain_fasta, 'fasta'):
        id = record.name
        value = len(record.seq)
        strain_dict[id] = value

    as_query_df['qcoverage'] = (as_query_df['length']/as_query_df['qseqid'].map(strain_dict))*100
    as_query_df['scoverage'] = (as_query_df['length']/as_query_df['sseqid'].map(MG1655_dict))*100

    as_target_df['qcoverage'] = (as_target_df['length']/as_target_df['qseqid'].map(MG1655_dict))*100
    as_target_df['scoverage'] = (as_target_df['length']/as_target_df['sseqid'].map(strain_dict))*100

    as_query_filtered = as_query_df.query("(evalue < 1e-5) & (pident >= 35) & (qcoverage >=70) & (scoverage >= 70)")
    as_target_filtered = as_target_df.query("(evalue < 1e-5) & (pident >= 35) & (qcoverage >=70) & (scoverage >= 70)")

    as_query_summary = as_query_filtered.groupby('qseqid').size().to_frame('num')
    as_target_summary = as_target_filtered.groupby('qseqid').size().to_frame('num')

    print('As query:')
    print(f"MAX: {as_query_summary['num'].max()}, MEAN: {as_query_summary['num'].mean()}, MED: {as_query_summary['num'].median()}, NUM>1: {len(as_query_summary.query('num >1'))}, NUM: {len(as_query_summary)}")

    print('As target:')
    print(f"MAX: {as_target_summary['num'].max()}, MEAN: {as_target_summary['num'].mean()}, MED: {as_target_summary['num'].median()}, NUM>1: {len(as_target_summary.query('num >1'))}, NUM: {len(as_target_summary)}")