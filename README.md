# Protisplice

A Python package for extracting splice site sequences from any organism given
a FASTA and GFF3 annotation file.

Protisplice allows for straightforward dataset generation for splice site classifiers.

In addition to locating true splice sites of any length, users can also generate
a "decoy" set of sequences that can be sampled from regions near, but not containing,
splice sites, or from intergenic regions. Users can define a buffer size, which
will start decoy sampling a defined number of bases away from the splice site.

Functions for calculating PPM, PWM, and PFM of the extracted sequences
are included in `protisplice.motif_scoring`.

Finally, functions for augmenting true splice sites and decoys can be found in
`protisplice.feature_eng`. These functions can mutate the consensus sequences
from true splice sites, as well as create consensus motifs for decoys. These
functions may help in creating more difficult negatives for model training,
which may improve the detection of false positives.

## Features

- Extract true and decoy splice site sequences from GFF3 and FASTA files of any organism
- Configurable extraction parameters (exon/intron lengths, buffer sizes)
- Output extracted sequences to FASTA
- Filtering transcripts on kallisto count data
- Functions for calculating position weight, frequency, and probability matrices
- Augmentation to strip or add consensus splice motifs

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

For a more in-depth tutorial, please see the included Jupyter notebook,
example_workflow.ipynb.

```python
from protisplice import SpliceSeqExtractor, ExtractionParams

# Basic usage
extractor = SpliceSeqExtractor(
    gff_path="annotations.gff3",
    fasta_path="genome.fasta"
)

true_ss = extractor.extract_splice_sites()

exon_decoys = extractor.sample_exonic_regions()
intronic_decoys = extractor.sample_intronic_regions()
intergenic_seqs = extractor.sample_intergenic_regions()
)

count = write_sequences(true_ss, "true_ss.fasta")
```

By default, `SpliceSeqExtractor` will find all available true splice sites, taking 40
bases from the exon and 80 bases from the intron for default sequence length of 120.

Additionally, the `SpliceSeqExtractor` class supports finding 'decoy' sequences by
sampling introns and exons from within a transcript, and will attempt to return a sequence of 
the same length as the true splice sites. 

Lastly, users can define a `buffer_size` which will prohibit 
sampling of introns/exons within `buffer_size` bases alongside either exon boundary. Note that when the `buffer_size` 
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

I was planning to add in some way of parsing VCFs in the future for detection of disrupted
splice sites, so expression filtering may be useful in combination with this if I decide
to go forward with that!

## File Format Requirements

### GFF3 File
- Standard GFF3 format with exon features
- Must contain transcript/mRNA entries with proper Parent relationships
- Biotype information recommended for filtering

### FASTA File
- Standard FASTA format
- If the file is not indexed, a .fai index will be created with pysam
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

For decoy sequences, the `junction_type` will be either
`intron`, `exon`, or `intergenic`.

## License

MIT License - see LICENSE file for details.

## Support

For questions or issues, please open an issue on the GitHub repository or contact bambha.m@northeastern.edu.