# Splice Site Prediction for Protist Genomes

## Description
This script identifies splice junctions (donor and acceptor sites) from a gene annotation file (GFF3 format) and extracts the flanking genomic sequences from a corresponding reference genome (FASTA format). It calculates a window around each splice junction coordinate, defined by a specified number of bases into the exon and intron regions.

Strand information is taken into account when outputting the sequences. For sequences on the (-) strand, the reverse complement will be returned. Additionally, intronic regions of the same window length are also extracted to serve as "decoy" samples.

This script was originally built as a data mining tool for downstream ML workflows.

## Features
- **Positive sample extraction**: True splice sites (donor and acceptor)
- **Negative sample extraction**: Random intronic regions as decoy sequences
- **Strand-aware processing**: Automatic reverse complement for negative strand sequences
- **Configurable windows**: Customizable exon and intron flanking regions

## Dependencies
* **Python 3.6 or higher**
* `pysam`: library for processing common bioinformatics data types
* `samtools`: for indexing FASTA files
* `biopython`: for FASTA I/O and reverse complement
* `pyyaml`: for loading YAML files containing test data (testing only)
* `pytest`: for running the tests (testing only)

## Installation

### Via pip
```bash
pip install pysam biopython pytest pyyaml
```

### Via conda (recommended)
```bash
conda install -c conda-forge -c bioconda pysam biopython samtools pytest pyyaml
```

### samtools installation
`samtools` is required for FASTA indexing:
```bash
# Via conda (recommended)
conda install -c conda-forge -c bioconda samtools

# Or follow instructions at: https://www.htslib.org/download/
```

**Important**: Ensure that the input GFF and FASTA files are from the same reference versions, otherwise junction locations may be incorrect.

## Usage

### Basic Usage
```bash
python3 extract_splice_seqs.py \
    --gff annotation.gff3 \
    --fasta genome.fa \
    --out1 positive_sequences.fa \
    --out2 negative_sequences.fa
```

### Full Command Line Options
```bash
python3 extract_splice_seqs.py \
    --gff <path-to-gff> \
    --fasta <path-to-fasta> \
    --out1 <positive-output> \
    --out2 <negative-output> \
    --n_exon <bases> \
    --n_intron <bases> \
    --buffer <buffer-size> \
    --verbose
```

### Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gff` | Required | Path to input GFF3 annotation file |
| `--fasta` | Required | Path to input FASTA genome file (must be indexed) |
| `--out1` | Required | Output file for positive (true splice site) sequences |
| `--out2` | Required | Output file for negative (decoy) sequences |
| `--n_exon` | 40 | Number of bases to include from exon region |
| `--n_intron` | 80 | Number of bases to include from intron region |
| `--buffer` | 50 | Buffer size around splice sites for negative sampling |
| `--verbose` | False | Enable detailed logging output |

## Input File Requirements

### GFF3 File Format
- Must contain `exon` features with proper transcript hierarchy
- Requires `Parent` or `transcript_id` attributes for transcript grouping
- Should follow standard GFF3 formatting conventions
- May add support for other types (GTF) later.

Example GFF3 entry:
```
chr1    source    exon    1000    1200    .    +    .    Parent=transcript:TRANSCRIPT_001
chr1    source    exon    1500    1700    .    +    .    Parent=transcript:TRANSCRIPT_001
```

### FASTA File Requirements
- Standard nucleotide FASTA format
- **Must be indexed** with `samtools faidx` before running
- Chromosome/contig names must match those in the GFF3 file

## Output Format

### Positive Sequences (True Splice Sites)
```
>chr1_donor_+_950_1120
ATCGATCGATCG...GTAAGT...TCGATCGATCG
>chr1_acceptor_+_1420_1580  
ATCGATCGATCG...TTGCAG...TCGATCGATCG
```

### Negative Sequences (Decoy/Intronic)
```
>chr1_intron_+_1250_1410
ATCGATCGATCGATCGATCGATCGATCGATCG...
```

Header format: `>{chromosome}_{type}_{strand}_{start}_{end}`

## Example Workflow

1. **Prepare your files:**
   ```bash
   # Index your FASTA file
   samtools faidx genome.fa
   
   # Verify GFF3 format
   head -20 annotation.gff3
   ```

2. **Run extraction:**
   ```bash
   python3 extract_splice_seqs.py \
       --gff annotation.gff3 \
       --fasta genome.fa \
       --out1 splice_sites_positive.fa \
       --out2 splice_sites_negative.fa \
       --n_exon 50 \
       --n_intron 100 \
       --verbose
   ```

3. **Check results:**
   ```bash
   # Count sequences
   grep -c ">" splice_sites_positive.fa
   grep -c ">" splice_sites_negative.fa
   
   # Verify sequence lengths
   head -4 splice_sites_positive.fa
   ```

## Troubleshooting

### Common Issues
1. **"FASTA index not found"**: Run `samtools faidx genome.fa`
2. **"No sequences extracted"**: Check GFF3 format and chromosome name matching
3. **Memory issues**: Use smaller buffer sizes or process in batches for very large genomes

## Testing
Run the test suite to verify installation:
```bash
pytest tests/
```

## Contact
- **Author**: Michael Bambha
- **Email**: bambha.m@northeastern.edu

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.