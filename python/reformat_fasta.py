from Bio import SeqIO
import pandas as pd
import re

genomes = ['Sq4a', 'GAOx1', 'MG1655']
for genome in genomes:
    input_fasta = f"../genomes/{genome}/cds/{genome}_cds.fna"
    input_protein = f"../genomes/{genome}/proteins/{genome}_proteins.faa"
    output_fasta = f"genome_fasta/{genome}_cds.fna"
    output_tsv = f"genome_fasta/{genome}_metadata.tsv"
    output_protein = f"genome_fasta/{genome}_proteins.faa"

    metadata = []
    new_records = []
    new_proteins = []

    for record in SeqIO.parse(input_fasta, "fasta"):

        desc = record.description

        attrs = dict(re.findall(r'\[([^=\]]+)=([^\]]+)\]', desc))

        locus_tag = attrs.get("locus_tag")

        if locus_tag is None:
            raise ValueError(f"Missing locus_tag in:\n{desc}")

        row = {"locus_tag": locus_tag}
        row.update(attrs)
        metadata.append(row)

        record.id = locus_tag
        record.name = locus_tag
        record.description = ""

        new_records.append(record)

    SeqIO.write(new_records, output_fasta, "fasta")

    metadata_df = pd.DataFrame(metadata)
    
    metadata_df.to_csv(
        output_tsv,
        sep='\t',
        index=False
    )
    protein_to_locus = dict(
        zip(metadata_df["protein_id"], metadata_df["locus_tag"])
    )
  
    for record in SeqIO.parse(input_protein, "fasta"):

        desc = record.description
        protein_id = desc.split()[0]
        locus_tag = protein_to_locus.get(protein_id)

        if locus_tag is None:
            raise ValueError(f"Protein ID not found in metadata: {protein_id}")

        record.id = locus_tag
        record.name = locus_tag
        record.description = ""

        new_proteins.append(record)

    SeqIO.write(new_proteins, output_protein, "fasta")

    print(f"{genome}: {len(new_records)} CDS, {len(new_proteins)} proteins")


    print(f"Processed {len(new_records)} nucleotide sequences")
    print(f"Processed {len(new_proteins)} protein sequences")
