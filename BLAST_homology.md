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

## Downloading Genome Sequences

The RefSeq genome annotation features (.gtf), genome sequences, genomic coding sequences, and proteins were downloaded from NCBI for Cablleronia str. [GAOx1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023631065.1/) and [Sq4a](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023170545.1/), and for E. coli str. K-12 substr. [MG1655](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/). These were transfered to the Emory Biology Server. 

To simplify the BLAST process, sequence identifiers need to be reformatted from their standard RefSeq formats. GEM iJO1366 uses locus-tags in their gene-reaction associations, so the protein fasta and cds fasta files must be reformatted so the locus tag is the identifier.

**Bash**
```
(microbial_GEMs) aenucko@bio:/data/ngerard/GEMs/microbial/BLAST$ python3 reformat_fasta.py
```
reformat_fasta.py can be found in the python folder of this repository. This code outputs to the directory aenucko@bio:/data/ngerard/GEMs/microbial/BLAST/genome_fasta/. For each strain (GAOx1, Sq4a, and MG1655), two fasta files are created (one containing cds and the other containing protein sequences, each identified with gene locus tags) and a metadata tsv file containing all of the gene information associated with the locus tag (gene, protein, protein_id, etc).

## Making BLAST Databases

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

## Running BLAST Searches

I wrote a bash script, named blast_search.sh, and ran it with SLURM to perform my BLAST searches. My bash script looked like this:

**Bash**
```
#!/bin/bash
#SBATCH --partition=day-long

source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbial_GEMs

blastp -query genome_fasta/Sq4a_proteins.faa -db blastdb/MG1655_blastdb/MG1655_blastdb -outfmt 6 -out blast_results/Sq4a/Sq4a_against_MG1655
blastp -query genome_fasta/MG1655_proteins.faa -db blastdb/Sq4a_blastdb/Sq4a_blastdb -outfmt 6 -out blast_results/Sq4a/MG1655_against_Sq4a

blastp -query genome_fasta/GAOx1_proteins.faa -db blastdb/MG1655_blastdb/MG1655_blastdb -outfmt 6 -out blast_results/GAOx1/GAOx1_against_MG1655
blastp -query genome_fasta/MG1655_proteins.faa -db blastdb/GAOx1_blastdb/GAOx1_blastdb -outfmt 6 -out blast_results/GAOx1/MG1655_against_GAOx1

conda deactivate
```
To run the script, I ran:

**Bash**
```
sbatch blast_search.sh
```
This code results in a table like the one below outputted to the location specified with the -out argument.

| qseqid | sseqid | pident | length | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| NCS66_RS00300 | b3251 | 71.143 | 350 | 95 | 3 | 1 | 347 | 1 | 347 | 0.0 | 506 |
| NCS66_RS00300 | b2526 | 27.126 | 247 | 139 | 9 | 15 | 223 | 24 | 267 | 1.65e-13 | 68.2 |
| NCS66_RS00300 | b0650 | 27.358 | 212 | 128 | 7 | 13 | 203 | 8 | 214 | 2.35e-11 | 61.6 |

## Summarizing BLAST hits

This table contains all of the BLAST hits for every query sequence. For model construction, I need to fine tune the level of scrutiny with which I will "accept" BLAST hits in order to minimize necessary gap filling while maintaining accuracy for what reactions the Caballeronia strains are capable of.

The paper referenced above uses the following criteria for BLAST hits: 
* At least 35% amino acid sequence identity
* E value less than 1e-5
* At least 70% coverage of both query and subject

Additionally, the paper uses reciprocal BLAST hits, meaning the match must be a hit with MG1655 as the subject and the query. 

For some basic statistics regarding BLAST hits, I wrote a python script that filters the BLAST hits for the above criteria, then counts how many hits each query sequence got (both with Cab and E. coli as the query). It then outputs the max number of hits for a query, the mean number of hits, the median number of hits, the number of query sequences with more than one hit, and the total number of query sequences with hits.

**Bash**
```
python3 count_hits.py
```
The output: 
> Sq4a: 
>> As query:
>>> MAX: 18, MEAN: 1.5998766954377313, MED: 1.0, NUM>1: 380, NUM: 1622
>>
>> As target:
>>> MAX: 19, MEAN: 1.777394900068918, MED: 1.0, NUM>1: 381, NUM: 1451
> 
> GAOx1:
>>As query:
>>> MAX: 17, MEAN: 1.6579945799457994, MED: 1.0, NUM>1: 453, NUM: 1845
>>
>> As target:
>>> MAX: 30, MEAN: 1.9921824104234527, MED: 1.0, NUM>1: 432, NUM: 1535

The script will also output the filtered tables as tsv files. 

## Constructing Draft Model from Reciprocal BLAST Hits

In order to use BLAST homology to pull reactions from the iJO1366 GEM, I need to create a list of MG1655 genes with reciprocal hits that were not filtered out, and determine which ones have corresponding reactions in iJO1366. Through manual inspection of that list, I will determine whether the filters need to be adjusted. 

match_hits.py is the python script I wrote. It outputs a tsv file containing all reciprocal BLAST hits and a tsv file containing the draft reconstruction reactions (GPRs not updated from iJO1366). Additionally, the number of reactions included for each draft reconstruction is printed. 

**Bash**
```
python3 match_hits.py
```

Output:
> Sq4a: 
>>total reactions kept: 891
>
>GAOx1: 
>>total reactions kept: 934

## Assigning Reactions to Orphan Genes

BLAST homology search assigned 891 total reactions for Sq4a and 934 total reactions for GAOx1. However, the majority of genes for each strain are not assigned reactions, called orphan genes. While many of these are structural or signal proteins, there are many which may catalyze reactions but have no direct homolog in MG1655. Firstly, using the python script, gtf_to_tsv.py, I created a tsv file containing only the orphan genes. 

**Bash**
```
(microbial_GEMs) aenucko@bio:/data/ngerard/GEMs/microbial/BLAST$ python3 gtf_to_tsv.py
```

Next, we need to go about assigning reactions to these orphan genes. Eggnog-mapper is an annotation software that assigns functional annotations to microbial genomes, including EC numbers, KEGG orthology terms, and KEGG reactions ids associated with the KO terms. Using the eggnog annotation csv file (eggnog annotations performed in Galaxy), for GAOx1 and Sq4a, I extracted KEGG reaction ids for orphan genes when they were present and mapped them to BiGG reaction ids, by cross referencing through MetaNetX. The script that performs this task is egg_to_rxns.py

For each strain, this script outputs three files: strain_reaction_set.txt containing a list of all BiGG reaction ids found from orphan genes, strain_gene_set.txt containing a list of all genes found to be associated with these reactions, and strain_genes_reactions.json containing a dictionary where each key is a BiGG reaction id, and each value is the list of genes associated with that reaction.

**Bash**
```
(microbial_GEMs) aenucko@bio:/data/ngerard/GEMs/microbial/orphan_genes$ python3 egg_to_rxns.py
```
Output:
> Sq4a:
>> NUMBER OF GENES: 284
>> 
>> NUMBER OF REACTIONS: 465
>
> GAOx1:
>> NUMBER OF GENES: 382
>> 
>> NUMBER OF REACTIONS: 544

## Updating Draft Reconstruction and Organizing Gap Genes




