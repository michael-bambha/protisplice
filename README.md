# Splice Site Prediction for Protist Genomes

## Description
This script identifies splice junctions (donor and acceptor sites) from a gene annotation file (GFF3 format) and extracts the flanking genomic sequences from a corresponding reference genome (FASTA format). It calculates a window around each splice junction coordinate, defined by a specified number of bases into the exon and intron regions.

Strand information is taken into account when outputting the sequences. For sequences on the (-) strand, the reverse complement will be returned. Additionally, intronic regions of the same window length are also extracted to serve as "decoy" samples.

The script supports filtering transcripts by biotype (e.g., protein-coding only) and expression levels, making it suitable for focused analysis of actively transcribed genes.

This script was originally built as a data mining tool for downstream ML workflows.

## Features
- **Positive sample extraction**: True splice sites (donor and acceptor)
- **Negative sample extraction**: Random intronic regions as decoy sequences
- **Strand-aware processing**: Automatic reverse complement for negative strand sequences
- **Configurable windows**: Customizable exon and intron flanking regions
- **Transcript filtering**: Filter by biotype (all transcripts or protein-coding only)
- **Expression-based filtering**: Filter transcripts by expression levels using RNA-seq quantification data
- **Multiple expression formats**: Support for Kallisto (Salmon and StringTie in future) output formats

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
    --transcript-filter <filter-type> \
    --expression-file <expression-file> \
    --min-expression <threshold> \
    --expression-format <format> \
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
| `--transcript-filter` | `protein_coding` | Include `all` or only `protein_coding` transcripts |
| `--expression-file` | None | Path to expression quantification file (TSV format) |
| `--min-expression` | None | Minimum expression threshold (TPM/FPKM) |
| `--expression-format` | `kallisto` | Format of expression file (`kallisto`, `salmon`, `stringtie`) |
| `--verbose` | False | Enable detailed logging output |

## Input File Requirements

### GFF3 File Format
- Must contain `exon` features with proper transcript hierarchy
- Must contain `mRNA` features with `biotype` attributes for transcript filtering
- Requires `Parent` or `transcript_id` attributes for transcript grouping
- Should follow standard GFF3 formatting conventions
- May add support for other types (GTF) later.

Example GFF3 entry:
```
chr1    source    mRNA     1000    2000    .    +    .    ID=transcript:TRANSCRIPT_001;biotype=protein_coding
chr1    source    exon     1000    1200    .    +    .    Parent=transcript:TRANSCRIPT_001
chr1    source    exon     1500    1700    .    +    .    Parent=transcript:TRANSCRIPT_001
```

### FASTA File Requirements
- Standard nucleotide FASTA format
- **Must be indexed** with `samtools faidx` before running
- Chromosome/contig names must match those in the GFF3 file

### Expression File Requirements (Optional)
The script supports expression quantification files from popular RNA-seq tools:

#### Kallisto Format
TSV file with columns: `target_id`, `length`, `eff_length`, `est_counts`, `tpm`
```
target_id       length  eff_length      est_counts      tpm
TRANSCRIPT_001  1500    1350.5         100.0           5.2
TRANSCRIPT_002  2000    1850.3         50.0            2.1
```

#### Salmon Format (planned)
Similar TSV format with transcript quantifications

#### StringTie Format (planned)
GTF-style output with FPKM values

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

## Example Workflows

### Basic Extraction (All Transcripts)
```bash
python3 extract_splice_seqs.py \
    --gff annotation.gff3 \
    --fasta genome.fa \
    --out1 splice_sites_positive.fa \
    --out2 splice_sites_negative.fa \
    --transcript-filter all \
    --verbose
```

### Protein-Coding Transcripts Only
```bash
python3 extract_splice_seqs.py \
    --gff annotation.gff3 \
    --fasta genome.fa \
    --out1 splice_sites_positive.fa \
    --out2 splice_sites_negative.fa \
    --transcript-filter protein_coding \
    --verbose
```

### Expression-Filtered Extraction
```bash
# First, prepare your files
samtools faidx genome.fa

# Extract sequences from highly expressed transcripts only
python3 extract_splice_seqs.py \
    --gff annotation.gff3 \
    --fasta genome.fa \
    --out1 splice_sites_positive_expressed.fa \
    --out2 splice_sites_negative_expressed.fa \
    --transcript-filter protein_coding \
    --expression-file kallisto_output.tsv \
    --min-expression 1.0 \
    --expression-format kallisto \
    --n_exon 50 \
    --n_intron 100 \
    --verbose
```

### Custom Window Sizes
```bash
python3 extract_splice_seqs.py \
    --gff annotation.gff3 \
    --fasta genome.fa \
    --out1 splice_sites_positive.fa \
    --out2 splice_sites_negative.fa \
    --n_exon 60 \
    --n_intron 120 \
    --buffer 75 \
    --verbose
```

## Advanced Usage

### Filtering Strategy
The script applies filters in the following order:
1. **Biotype filtering**: Include only specified transcript types
2. **Expression filtering**: Include only transcripts above expression threshold
3. **Structural filtering**: Require transcripts with ≥2 exons for splice junction analysis

### Expression Threshold Selection
- **Low expression (0.1-1.0 TPM)**: Include more transcripts, may include noise
- **Moderate expression (1.0-5.0 TPM)**: Balance between coverage and quality
- **High expression (>5.0 TPM)**: Focus on highly active genes only

### Window Size Considerations
- **Smaller windows (20-40 nt)**: Focus on immediate splice site motifs
- **Larger windows (80-200 nt)**: Capture broader regulatory context
- **Buffer size**: Prevent negative samples from overlapping true splice sites

## Troubleshooting

### Common Issues
1. **"FASTA index not found"**: Run `samtools faidx genome.fa`
2. **"No sequences extracted"**: Check GFF3 format and chromosome name matching
3. **"Expression file not found"**: Verify path and file format
4. **Memory issues**: Use smaller buffer sizes or process in batches for very large genomes
5. **Transcript ID mismatch**: Check that transcript IDs in GFF3 match expression file

### Expression File Issues
- **Transcript ID formatting**: Script handles common variations (`transcript:ID`, version numbers)
- **Missing expression data**: Transcripts without expression data are assigned 0.0 TPM
- **Format detection**: Verify column headers match expected format

### Debugging Tips
```bash
# Check transcript extraction
python3 extract_splice_seqs.py --gff annotation.gff3 --fasta genome.fa \
    --out1 test_pos.fa --out2 test_neg.fa --verbose

# Verify expression file parsing
head -5 your_expression_file.tsv

# Check GFF3 transcript structure
grep -E "(mRNA|exon)" annotation.gff3 | head -10
```

## Testing
Run the test suite to verify installation:
```bash
pytest tests/
```

## Performance Notes
- Processing time scales with genome size and transcript count
- Memory usage is generally modest due to streaming processing
- Expression filtering can significantly reduce runtime for large annotations
- Will add batch processing of annotations + multithread support in the future 

## Contact
- **Author**: Michael Bambha
- **Email**: bambha.m@northeastern.edu

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.