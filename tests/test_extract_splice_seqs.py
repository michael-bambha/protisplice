"""
Test suite for extract_splice_seqs.py
"""

import os
from pathlib import Path
import tempfile
import pysam
import pytest
from load_test_cases import (
    load_grpexons_cases,
    load_gtf_test_cases_yaml,
    load_splicejxn_cases,
    load_intron_coords_cases,
    load_window_coords_cases,
    load_extract_seq_cases,
)
from extract_splice_seqs import (
    _parse_transcript_id,
    group_exons_by_transcript,
    get_splice_junctions,
    extract_sequence,
    get_window_coords,
    get_intron_coords,
)

TEST_CASES_YAML_PATH = Path("tests/test_data").resolve()


@pytest.mark.parametrize(
    "gtf_str, expected",
    load_gtf_test_cases_yaml(
        path=TEST_CASES_YAML_PATH / "test__parse_transcript_id.yaml"
    ),
)
def test__parse_transcript_id(gtf_str, expected):
    assert _parse_transcript_id(gtf_str) == expected


@pytest.mark.parametrize(
    "gtf_content, expected",
    load_grpexons_cases(
        path=TEST_CASES_YAML_PATH / "test_group_exons_by_transcript_cases.yaml"
    ),
)
def test_group_exons_by_transcript(gtf_content, expected):
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, encoding="utf-8", suffix=".gtf"
    ) as tmp_gtf:
        tmp_gtf.write(gtf_content)
        gtf_path = tmp_gtf.name
    try:
        result = group_exons_by_transcript(gtf_path)
        assert result == expected
    finally:
        os.remove(gtf_path)


@pytest.mark.parametrize(
    "transcripts, expected",
    load_splicejxn_cases(
        path=TEST_CASES_YAML_PATH / "test_get_splice_junctions_cases.yaml"
    ),
)
def test_get_splice_junctions(transcripts, expected):
    result = get_splice_junctions(transcripts)
    result_sorted = sorted(result, key=lambda x: x["id"])
    expected_sorted = sorted(expected, key=lambda x: x["id"])
    assert result_sorted == expected_sorted


@pytest.mark.parametrize(
    "transcripts, expected",
    load_intron_coords_cases(
        path=TEST_CASES_YAML_PATH / "test_get_intron_coords_cases.yaml"
    ),
)
def test_get_intron_coords(transcripts, expected):
    result = get_intron_coords(transcripts)
    assert result == expected


@pytest.mark.parametrize(
    "strand_arg, junc_type_arg, coord_arg, n_exon_arg, n_intron_arg, expected_val",
    load_window_coords_cases(
        path=TEST_CASES_YAML_PATH / "test_get_window_coords_cases.yaml"
    ),
)
def test_get_window_coords(
    strand_arg, junc_type_arg, coord_arg, n_exon_arg, n_intron_arg, expected_val
):
    result = get_window_coords(
        strand_arg, junc_type_arg, coord_arg, n_exon_arg, n_intron_arg
    )
    assert result == expected_val


@pytest.mark.parametrize(
    "seq_id_arg, win_start_arg, win_end_arg, fasta_content_str, expected_seq_val",
    load_extract_seq_cases(path=TEST_CASES_YAML_PATH / "test_extract_seq_cases.yaml"),
)
def test_extract_sequence(
    seq_id_arg: str,
    win_start_arg: int,
    win_end_arg: int,
    fasta_content_str: str,
    expected_seq_val: str,
    tmp_path: Path,
):
    """
    Tests the extract_sequence function using test cases loaded from YAML.

    Args:
        seq_id_arg (str): Sequence ID for the desired sequence.
        win_start_arg (int): Start coordinate (1-based, inclusive).
        win_end_arg (int): End coordinate (1-based, inclusive).
        fasta_content_str (str): The string content of the FASTA file for this test case.
        expected_seq_val (str): The expected sequence string.
        tmp_path (Path): Pytest fixture for a temporary directory.
    """
    fasta_file_path = tmp_path / "test.fa"
    actual_sequence = ""
    fasta_file_path.write_text(fasta_content_str)
    if not fasta_content_str.strip() or not any(
        line.startswith(">") for line in fasta_content_str.splitlines()
    ):
        if expected_seq_val == "":
            if fasta_file_path.exists() and fasta_file_path.read_text().strip():
                pysam.faidx(str(fasta_file_path))
            with pysam.FastaFile(str(fasta_file_path)) as fasta_obj:
                actual_sequence = extract_sequence(
                    fasta_obj, seq_id_arg, win_start_arg, win_end_arg
                )
        assert actual_sequence == expected_seq_val
