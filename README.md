# Splice Site Prediction for Protist Genomes

## Description

This script identifies splice junctions (donor and acceptor sites) from a gene annotation file (GTF format) and extracts the flanking genomic sequences from a corresponding reference genome (FASTA format). It calculates a window around each splice junction coordinate, defined by a specified number of bases into the exon and intron regions.

Strand information is taken into account when outputting the sequences. For sequences on the (-) strand,
the reverse complement will be returned.

This script was originally built as a data mining tool for downstream ML workflows. 

## Dependencies

* python 3.6 or higher
* `pysam`: library for processing common bioinformatics data types.
* `samtools`: for indexing FASTA files
* `pyyaml`: for loading YAML files containing test data
* `pytest`: for running the tests

## Installation

pysam can be installed via pip:

```bash
pip install pysam
```

A GTF annotation file and both a FASTA and indexed FASTA file are required to execute the script.

To index a fasta file:

```bash
samtools faidx yourfile.fa
```

## Usage

```bash
python3 splice_sites_protists.py --gtf <path-to-gtf> --fasta <path-to-fasta> --n_exon <bases> --n_intron <bases>
```

