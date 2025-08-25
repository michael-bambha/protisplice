# Protisplice

A Python package for extracting splice site sequences from any organism given
a FASTA and GFF3 annotation file.

Protisplice allows for straightforward dataset generation for splice site classifiers.
In addition to locating true splice sites of any length, users can also generate
a "decoy" set of sequences that can be sampled from regions near, but not containing,
splice sites, or from intergenic regions. Users can define a buffer size, which
will start decoy sampling a defined number of bases away from the splice site.

Lastly, functions for calculating PPM, PWM, and PFM of the extracted sequences
are included in `protisplice.motif_scoring`.

## Features

- Extract true and decoy splice site sequences from GFF3 and FASTA files of any organism
- Configurable extraction parameters (exon/intron lengths, buffer sizes)
- Output extracted sequences to FASTA
- Filtering transcripts on kallisto count data
- Functions for calculating position weight, frequency, and probability matrices

## Installation

### From Source

```bash
git clone https://github.com/michael-bambha/protisplice.git
cd protisplice
pip install -e .
```

### Requirements

- Python ≥ 3.8
- pysam ≥ 0.19.0
- biopython ≥ 1.79
- pandas ≥ 1.4
- numpy ≥ 2.0
- pytest (dev)

## Quick Start

For an additional tutorial, please see the included Jupyter notebook,
example_workflow.ipynb.

```python
from protisplice import SpliceSeqExtractor, ExtractionParams

# Basic usage
extractor = SpliceSeqExtractor(
    gff_path="annotations.gff3",
    fasta_path="genome.fasta"
)

# Write to FASTA files
pos_count, neg_count = extractor.write_sequences_to_fasta(
    results,
    "positive_sequences.fasta",
    "negative_sequences.fasta"
)

print(f"Extracted {pos_count} positive and {neg_count} negative sequences")
```

This method will split true splice sites and decoy splice sites into two separate files. 
By default, `SpliceSeqExtractor` will find all available true splice sites, taking 40
bases from the exon and 80 bases from the intron for default sequence length of 120.

Additionally, the `SpliceSeqExtractor` class supports finding 'decoy' sequences by
sampling introns and exons from within a transcript, and will attempt to return a sequence of 
the same length as the true splice sites. Lastly, users can define a `buffer_size` which will prohibit 
sampling of introns within `buffer_size` bases alongside either exon boundary. Note that when the `buffer_size` 
is applied, certain introns or exons may too short to reach the desired window size, and will be discarded. 


## Advanced Usage

### Custom Parameters

```python
from protisplice import ExtractionParams, TranscriptFilter

params = ExtractionParams(
    n_exon=50,      # Bases from exon region
    n_intron=100,   # Bases from intron region  
    buffer_size=75  # Buffer around splice sites for decoy sampling
)

extractor = SpliceSeqExtractor(
    gff_path="annotations.gff3",
    fasta_path="genome.fasta",
    params=params,
    transcript_filter=TranscriptFilter.PROTEIN_CODING
)
```
Additionally, the `transcript_filter` argument can be passed in to filter out
certain transcripts. Currently, supported filters are `TranscriptFilter.ALL` and
`TranscriptFilter.PROTEIN_CODING`. By default, all transcripts are included,
but users can filter for only protein-coding transcripts as dictated by the GFF3
annotations.

### Expression Filtering

```python
extractor = SpliceSeqExtractor(
    gff_path="annotations.gff3",
    fasta_path="genome.fasta",
    expression_file="abundance.tsv",
    min_expression=1.0,
    expression_format="kallisto"
)
```
Lastly, transcripts can be filtered out by count data. Currently, only kallisto
format is supported for this, but more formats like Salmon will be added in the future.
It is strongly recommended that the expression data be normalized for most workflows.

### Step-by-Step Extraction

```python
# Extract positive and negative sequences separately
true_ss  = extractor.extract_positive_sequences()
decoy_ss = extractor.extract_negative_sequences(len(positive_sequences))

# Access raw data
transcripts = extractor.transcripts
junctions = extractor.junctions
info = extractor.get_info()
print(info)
```

## File Format Requirements

### GFF3 File
- Standard GFF3 format with exon features
- Must contain transcript/mRNA entries with proper Parent relationships
- Biotype information recommended for filtering

### FASTA File
- Standard FASTA format
- **Must be indexed** with `samtools faidx` before use
- Sequence IDs must match those in GFF3 file

### Expression Files

**Kallisto format:**
```
target_id	length	eff_length	est_counts	tpm
transcript1	1000	800	150.5	12.3
```

**GTEx format:**
```
# GTEx expression matrix format
gene_id	gene_name	tissue1	tissue2	tissue3
ENSG00000001	GENE1	5.2	3.1	8.7
```

## Output Format

FASTA headers contain sequence metadata:
```
>chr1_donor_+_12345_12465
ATCGATCGATCG...

>chr2_acceptor_-_67890_68010
GCTAGCTAGCTA...
```

Header format: `>{seqid}_{junction_type}_{strand}_{start}_{end}`


## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Support

For questions or issues, please open an issue on the GitHub repository or contact bambha.m@northeastern.edu.