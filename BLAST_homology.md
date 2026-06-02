# Create Genome-Scale Metabolic Model Using BLAST Homology Search

## Description

Follows microbial metabolic-model reconstruction protocol outlined in [Ankrah, N. Y. D., Luan, J., and Douglas, A. E. (2017)](https://doi.org/10.1128/jb.00872-16).

I performed reciprocal BLAST searches of the Caballeronia genomes against E. coli str. K-12 substr. MG1655, determining orthologous genes. In the case of bidirectional gene matches, the associated reactions from the E. coli str. K-12 substr. MG1655 GEM, iJO1366, were included in the Caballeronia metabolic-model reconstruction if: E value < 1e-5, AA sequence identity > 35%, and match length > 70% for both subject and query sequences.

This process was completed on the Emory Biology Server.

## Installing BLAST+ Suite

First, create a conda environment.

**Bash**
```
conda create -n microbial_GEMs
conda activate microbial_GEMs
```
Then, install blast and Bio from bioconda.

**Bash**
```
conda install -c bioconda blast
conda install -c bioconda Bio
```

## Downloading Genome Sequences and Making BLAST Databases

The RefSeq genome annotation features (.gtf), genome sequences, genomic coding sequences, and proteins were downloaded from NCBI for Cablleronia str. [GAOx1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023631065.1/) and [Sq4a](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023170545.1/), and for E. coli str. K-12 substr. [MG1655](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/). These were transfered to the Emory Biology Server. 

To simplify the BLAST process, sequence identifiers need to be reformatted from their standard RefSeq formats. GEM iJO1366 uses locus-tags in their gene-reaction associations, so the protein fasta and cds fasta files must be reformatted so the locus tag is the identifier.

**Bash**
```
(microbial_GEMs) aenucko@bio:/data/ngerard/GEMs/microbial/BLAST$ python3 reformat_fasta.py
```
reformat_fasta.py can be found in the python folder of this repository. This code outputs to the directory aenucko@bio:/data/ngerard/GEMs/microbial/BLAST/genome_fasta/. For each strain (GAOx1, Sq4a, and MG1655), two fasta files are created (one containing cds and the other containing protein sequences, each identified with gene locus tags) and a metadata tsv file containing all of the gene information associated with the locus tag (gene, protein, protein_id, etc).

To make the BLAST databases for each microbial strain, run the following commands:

**Bash**
```{bash}
(microbial_GEMs) aenucko@bio:/data/ngerard/GEMs/microbial/BLAST$ makeblastdb \
    -in genome_fasta/GAOx1_proteins.faa \
    -title "Caballeronia sp. GAOx1 BLAST database" \
    -out blastdb/GAOx1_blastdb/GAOx1_blastdb \
    -dbtype prot

(microbial_GEMs) aenucko@bio:/data/ngerard/GEMs/microbial/BLAST$ makeblastdb \
    -in genome_fasta/Sq4a_proteins.faa \
    -title "Caballeronia sp. Sq4a BLAST database" \
    -out blastdb/Sq4a_blastdb/Sq4a_blastdb \
    -dbtype prot

(microbial_GEMs) aenucko@bio:/data/ngerard/GEMs/microbial/BLAST$ makeblastdb \
    -in genome_fasta/MG1655_proteins.faa \
    -title "E. coli str. K-12 substr. MG1655 BLAST database" \
    -out blastdb/MG1655_blastdb/MG1655_blastdb \
    -dbtype prot
```

Confirm the databases were created

**Bash**
```
cd blastdb/GAOx1_blastdb
ls
```
The output should look something like this:
> GAOx1_blastdb.pdb  GAOx1_blastdb.phr  GAOx1_blastdb.pin  GAOx1_blastdb.pot  GAOx1_blastdb.psq  GAOx1_blastdb.ptf  GAOx1_blastdb.pto









