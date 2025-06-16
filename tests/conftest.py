"""
Shared pytest fixtures and configuration for splice sequence extractor tests.
"""
# pylint:disable=redefined-outer-name

import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock
import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


@pytest.fixture(scope="session")
def test_data_dir():
    """Create a temporary directory for test data that persists for the session"""
    temp_dir = tempfile.mkdtemp(prefix="splice_test_")
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing"""
    return """>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"""


@pytest.fixture
def create_test_files(test_data_dir, sample_fasta_content):
    """Create actual test files with content"""

    def _create_files(gff_content="", fasta_content=None):
        if fasta_content is None:
            fasta_content = sample_fasta_content

        # Create files
        gff_file = test_data_dir / "test.gff3"
        fasta_file = test_data_dir / "test.fa"
        fai_file = test_data_dir / "test.fa.fai"

        # Write content
        with open(gff_file, "w", encoding='utf-8') as f:
            f.write(gff_content)

        with open(fasta_file, "w", encoding='utf-8') as f:
            f.write(fasta_content)

        # Create simple fai file
        with open(fai_file, "w", encoding='utf-8') as f:
            f.write("chr1\t1000\t6\t70\t71\n")
            f.write("chr2\t1000\t1006\t70\t71\n")

        return {"gff": str(gff_file), "fasta": str(fasta_file), "fai": str(fai_file)}

    return _create_files


@pytest.fixture
def mock_logger():
    """Mock logger for testing"""
    return Mock()
