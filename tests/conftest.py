"""
Shared pytest fixtures and configuration for splice sequence extractor tests.
"""

# pylint:disable=redefined-outer-name

import tempfile
from pathlib import Path
import pytest
import pysam

from protisplice.data_models import (
    ExtractionParams,
    SpliceJunction,
    JunctionType,
    StrandType,
    Transcript,
    TranscriptInfo
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
def sample_transcript():
    """Sample transcript for testing"""
    return Transcript(
        info=TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE),
        exons=[(100, 200), (300, 400), (500, 600)],
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
chr1\tensembl\tmRNA\t100\t600\t.\t+\t.\tID=transcript1;biotype=protein_coding
chr1\tensembl\texon\t100\t200\t.\t+\t.\tParent=transcript1
chr1\tensembl\texon\t300\t400\t.\t+\t.\tParent=transcript1
chr1\tensembl\texon\t500\t600\t.\t+\t.\tParent=transcript1
chr2\tensembl\tmRNA\t100\t300\t.\t-\t.\tID=transcript2;biotype=protein_coding
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
