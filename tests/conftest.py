"""
Shared pytest fixtures and configuration for splice sequence extractor tests.
Updated to include new data models and functionality.
"""

# pylint:disable=redefined-outer-name

import tempfile
from pathlib import Path
import pytest
import pysam

from protisplice.data_models import (
    ExtractionParams,
    SamplingParams,
    SpliceJunction,
    JunctionData,
    JunctionType,
    StrandType,
    Transcript,
    TranscriptInfo,
    Gene,
)


@pytest.fixture
def temp_dir():
    """Create temporary directory"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def extraction_params():
    """Extraction parameters for testing"""
    return ExtractionParams(n_exon=40, n_intron=80, buffer_size=50)


@pytest.fixture
def sampling_params():
    """Sampling parameters for testing"""
    return SamplingParams(window_size=120, buffer_size=50)


@pytest.fixture
def sample_transcript():
    """Sample transcript for testing"""
    return Transcript(
        info=TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE),
        exons=[(100, 200), (300, 400), (500, 600)],
    )


@pytest.fixture
def sample_negative_transcript():
    """Sample transcript on negative strand for testing"""
    return Transcript(
        info=TranscriptInfo(seqid="chr2", strand=StrandType.NEGATIVE),
        exons=[(100, 150), (250, 300)],
    )


@pytest.fixture
def sample_junctions():
    """Sample splice junctions for testing"""
    return [
        SpliceJunction(
            id="transcript1_donor_0",
            seqid="chr1",
            coord=200,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        ),
        SpliceJunction(
            id="transcript1_acceptor_1",
            seqid="chr1",
            coord=300,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.ACCEPTOR,
        ),
    ]


@pytest.fixture
def sample_genes():
    """Sample genes for testing"""
    return {
        "gene1": Gene(
            gene_id="gene1",
            seq_id="chr1",
            start=100,
            end=600,
            strand=StrandType.POSITIVE,
        ),
        "gene2": Gene(
            gene_id="gene2",
            seq_id="chr2",
            start=100,
            end=300,
            strand=StrandType.NEGATIVE,
        ),
    }


@pytest.fixture
def sample_junction_data():
    """Sample junction data for testing"""
    junction = SpliceJunction(
        id="test_junction",
        seqid="chr1",
        coord=200,
        strand=StrandType.POSITIVE,
        junction_type=JunctionType.DONOR,
    )
    return JunctionData(
        junction=junction,
        window_start=160,
        window_end=280,
        sequence="A" * 40 + "GT" + "T" * 78,  # 40 exon + GT + 78 intron = 120bp
    )


@pytest.fixture
def sample_fasta():
    """Sample FASTA for testing"""
    return {
        "chr1": "A" * 100
        + "GTAAGT"
        + "T" * 194
        + "TTGCAG"
        + "C" * 194
        + "GTAAGT"
        + "G" * 194
        + "TTGCAG"
        + "A" * 100,
        "chr2": "C" * 500,
    }


@pytest.fixture
def test_fasta_file(temp_dir, sample_fasta):
    """Create test FASTA file"""
    fasta_path = temp_dir / "test.fasta"
    with open(fasta_path, "w", encoding="utf-8") as f:
        for seq_id, sequence in sample_fasta.items():
            f.write(f">{seq_id}\n{sequence}\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


@pytest.fixture
def test_gff_content():
    """Sample GFF3 content for testing"""
    return """##gff-version 3
chr1\tensembl\tgene\t100\t600\t.\t+\t.\tID=gene:gene1;biotype=protein_coding
chr1\tensembl\tmRNA\t100\t600\t.\t+\t.\tID=transcript1;Parent=gene:gene1;biotype=protein_coding
chr1\tensembl\texon\t100\t200\t.\t+\t.\tParent=transcript1
chr1\tensembl\texon\t300\t400\t.\t+\t.\tParent=transcript1
chr1\tensembl\texon\t500\t600\t.\t+\t.\tParent=transcript1
chr2\tensembl\tgene\t100\t300\t.\t-\t.\tID=gene:gene2;biotype=protein_coding
chr2\tensembl\tmRNA\t100\t300\t.\t-\t.\tID=transcript2;Parent=gene:gene2;biotype=protein_coding
chr2\tensembl\texon\t100\t150\t.\t-\t.\tParent=transcript2
chr2\tensembl\texon\t250\t300\t.\t-\t.\tParent=transcript2
"""


@pytest.fixture
def test_gff_file(temp_dir, test_gff_content):
    """Create a test GFF3 file"""
    gff_path = temp_dir / "test.gff3"
    with open(gff_path, "w", encoding="utf-8") as f:
        f.write(test_gff_content)
    return gff_path


@pytest.fixture
def test_expression_data():
    """Sample expression data for testing"""
    return {"transcript1": 5.0, "transcript2": 0.5, "transcript3": 2.0}


@pytest.fixture
def test_kallisto_file(temp_dir):
    """Create a test Kallisto expression file"""
    kallisto_path = temp_dir / "abundance.tsv"
    content = """target_id\tlength\teff_length\test_counts\ttpm
transcript1\t1000\t800\t100.0\t5.0
transcript2\t800\t600\t10.0\t0.5
transcript3\t1200\t1000\t50.0\t2.0
"""
    with open(kallisto_path, "w", encoding="utf-8") as f:
        f.write(content)
    return kallisto_path


@pytest.fixture
def test_gtex_file(temp_dir):
    """Create a test GTEx expression file"""
    gtex_path = temp_dir / "gtex_expression.txt"
    content = """# GTEx Expression Data
# Version: V8
# Sample metadata
gene_id\tgene_name\tsample1\tsample2\tsample3
gene1\tGENE1\t5.2\t3.1\t8.7
gene2\tGENE2\t0.5\t0.3\t0.8
gene3\tGENE3\t2.1\t4.2\t1.8
"""
    with open(gtex_path, "w", encoding="utf-8") as f:
        f.write(content)
    return gtex_path


@pytest.fixture
def long_sequence():
    """Long sequence for comprehensive testing"""
    # Create a 1000bp sequence with known motifs
    sequence = (
        "A" * 200  # Exon 1
        + "GT"
        + "T" * 98
        + "AG"  # Intron 1 (100bp)
        + "C" * 200  # Exon 2
        + "GT"
        + "G" * 98
        + "AG"  # Intron 2 (100bp)
        + "T" * 200  # Exon 3
    )
    return sequence


@pytest.fixture
def complex_gff_content():
    """More complex GFF3 content for integration testing"""
    return """##gff-version 3
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:ENSG001;gene_name=TEST1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=ENST001.1;Parent=gene:ENSG001;transcript_biotype=protein_coding
chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tParent=ENST001.1
chr1\tensembl\texon\t2000\t2300\t.\t+\t.\tParent=ENST001.1
chr1\tensembl\texon\t3000\t3400\t.\t+\t.\tParent=ENST001.1
chr1\tensembl\texon\t4500\t5000\t.\t+\t.\tParent=ENST001.1
chr1\tensembl\tgene\t6000\t8000\t.\t-\t.\tID=gene:ENSG002;gene_name=TEST2;biotype=protein_coding
chr1\tensembl\tmRNA\t6000\t8000\t.\t-\t.\tID=ENST002.2;Parent=gene:ENSG002;transcript_biotype=protein_coding
chr1\tensembl\texon\t6000\t6500\t.\t-\t.\tParent=ENST002.2
chr1\tensembl\texon\t7000\t7300\t.\t-\t.\tParent=ENST002.2
chr1\tensembl\texon\t7800\t8000\t.\t-\t.\tParent=ENST002.2
chr2\tensembl\tgene\t1000\t3000\t.\t+\t.\tID=gene:ENSG003;gene_name=TEST3;biotype=lncRNA
chr2\tensembl\tlncRNA\t1000\t3000\t.\t+\t.\tID=ENST003.1;Parent=gene:ENSG003;transcript_biotype=lncRNA
chr2\tensembl\texon\t1000\t1500\t.\t+\t.\tParent=ENST003.1
chr2\tensembl\texon\t2500\t3000\t.\t+\t.\tParent=ENST003.1
"""


@pytest.fixture
def complex_gff_file(temp_dir, complex_gff_content):
    """Create a more complex test GFF3 file"""
    gff_path = temp_dir / "complex_test.gff3"
    with open(gff_path, "w", encoding="utf-8") as f:
        f.write(complex_gff_content)
    return gff_path


@pytest.fixture
def complex_fasta_file(temp_dir):
    """Create a complex test FASTA file with realistic sequences"""
    fasta_path = temp_dir / "complex_test.fasta"

    # Create sequences with known splice sites
    chr1_seq = (
        "N" * 1000  # Start padding
        + "A" * 200  # Exon 1 (1000-1200)
        + "GT"
        + "T" * 797
        + "AG"  # Intron 1 (1201-1999)
        + "C" * 300  # Exon 2 (2000-2300)
        + "GT"
        + "G" * 697
        + "AG"  # Intron 2 (2301-2999)
        + "T" * 400  # Exon 3 (3000-3400)
        + "GT"
        + "A" * 1098
        + "AG"  # Intron 3 (3401-4499)
        + "G" * 500  # Exon 4 (4500-5000)
        + "N" * 1000  # Gap
        + "C" * 500  # Gene 2 Exon 1 (6000-6500)
        + "GT"
        + "T" * 498
        + "AG"  # Gene 2 Intron 1 (6501-6999)
        + "A" * 300  # Gene 2 Exon 2 (7000-7300)
        + "GT"
        + "G" * 498
        + "AG"  # Gene 2 Intron 2 (7301-7799)
        + "T" * 200  # Gene 2 Exon 3 (7800-8000)
        + "N" * 1000  # End padding
    )

    chr2_seq = (
        "N" * 1000  # Start padding
        + "A" * 500  # Exon 1 (1000-1500)
        + "GT"
        + "T" * 998
        + "AG"  # Intron (1501-2499)
        + "C" * 500  # Exon 2 (2500-3000)
        + "N" * 1000  # End padding
    )

    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(f">chr1\n{chr1_seq}\n")
        f.write(f">chr2\n{chr2_seq}\n")

    pysam.faidx(str(fasta_path))
    return fasta_path
