"""
Functions for loading in test data from YAML
"""

import pytest
import yaml


def load_gtf_test_cases_yaml(path):
    """Loads test cases from a YAML file."""
    test_cases = []
    with open(path, "r", encoding="utf-8") as f:
        all_data = yaml.safe_load(f)
        for data in all_data:
            test_id = data.get("id", data["input_string"][:30])
            test_cases.append(
                pytest.param(
                    data["input_string"], data["expected_attributes"], id=test_id
                )
            )
    return test_cases


def load_grpexons_cases(path):
    with open(path, "r", encoding="utf-8") as f:
        cases = yaml.safe_load(f)
    test_cases = []
    for case in cases:
        test_id = case["test_id"]
        gtf_content = case.get("gtf_content", "")
        expected_output = case["expected_output"]
        for _, data in expected_output.items():
            if "exons" in data and isinstance(
                data["exons"], list
            ):  # convert list of lists to list of tuples
                data["exons"] = [tuple(exon) for exon in data["exons"]]
        test_cases.append(pytest.param(gtf_content, expected_output, id=test_id))
    return test_cases


def load_splicejxn_cases(path):
    with open(path, "r", encoding="utf-8") as f:
        cases = yaml.safe_load(f)
    test_cases = []
    for case in cases:
        test_id = case["test_id"]
        input_transcripts_raw = case.get("input_transcripts", {})
        expected_junctions = case.get("expected_junctions", [])
        input_transcripts_processed = {}
        for transcript_id, data in input_transcripts_raw.items():
            processed_exons = []
            if "exons" in data and isinstance(data["exons"], list):
                for exon_coords in data["exons"]:
                    if isinstance(exon_coords, list) and len(exon_coords) == 2:
                        processed_exons.append(tuple(exon_coords))
            input_transcripts_processed[transcript_id] = {
                "info": data.get("info", {}),
                "exons": processed_exons,
            }
        test_cases.append(
            pytest.param(input_transcripts_processed, expected_junctions, id=test_id)
        )
    return test_cases


def load_intron_coords_cases(path):
    with open(path, "r", encoding="utf-8") as f:
        cases = yaml.safe_load(f)
    test_cases = []

    for case in cases:
        test_id = case.get("test_id", "unnamed_intron_test")
        input_transcripts_raw = case.get("input_transcripts", {})
        expected_output_raw = case.get("expected_output", {})
        input_transcripts_processed = {}

        for transcript_id, data in input_transcripts_raw.items():
            processed_exons_input = []
            if "exons" in data and isinstance(data["exons"], list):
                for exon_coords in data["exons"]:
                    if isinstance(exon_coords, list) and len(exon_coords) == 2:
                        processed_exons_input.append(tuple(exon_coords))
            input_transcripts_processed[transcript_id] = {
                "info": data.get("info", {}),
                "exons": processed_exons_input,
            }

        expected_output_processed = {}
        for transcript_id, data in expected_output_raw.items():
            processed_exons_expected = []
            if "exons" in data and isinstance(data["exons"], list):
                for exon_coords in data["exons"]:
                    if isinstance(exon_coords, list) and len(exon_coords) == 2:
                        processed_exons_expected.append(tuple(exon_coords))
            processed_introns_expected = []
            if "introns" in data and isinstance(data["introns"], list):
                for intron_coords in data["introns"]:
                    if isinstance(intron_coords, list) and len(intron_coords) == 2:
                        processed_introns_expected.append(tuple(intron_coords))

            expected_output_processed[transcript_id] = {
                "info": data.get("info", {}),
                "exons": processed_exons_expected,
                "introns": processed_introns_expected,
            }

        test_cases.append(
            pytest.param(
                input_transcripts_processed, expected_output_processed, id=test_id
            )
        )
    return test_cases


def load_window_coords_cases(path):
    """
    Loads window coordinate test cases from a YAML file and prepares them for pytest.

    Args:
        yaml_path (Path): The path to the YAML file containing the test cases.

    Returns:
        List[Any]: A list of pytest.param objects, each containing the
                   individual input arguments (strand, junc_type, coord, n_exon, n_intron)
                   and the expected_window tuple, plus a test_id.
                   Returns an empty list if the file is not found or is empty.
    """
    with open(path, "r", encoding="utf-8") as f:
        raw_test_cases = yaml.safe_load(f)
    test_cases = []
    for case in raw_test_cases:
        test_id = case.get("test_id", "unnamed_window_test")
        inputs = case.get("inputs", {})
        strand = inputs.get("strand")
        junc_type = inputs.get("junc_type")
        coord = inputs.get("coord")
        n_exon = inputs.get("n_exon")
        n_intron = inputs.get("n_intron")
        expected_window_raw = case.get("expected_window")
        expected_window_tuple = None  # null is None when loaded w/ PyYAML
        if isinstance(expected_window_raw, list) and len(expected_window_raw) == 2:
            val1 = expected_window_raw[0]
            val2 = expected_window_raw[1]
            expected_window_tuple = (val1, val2)
        elif expected_window_raw is None:
            expected_window_tuple = (None, None)

        test_cases.append(
            pytest.param(
                strand,
                junc_type,
                coord,
                n_exon,
                n_intron,
                expected_window_tuple,
                id=test_id,
            )
        )
    return test_cases


def load_extract_seq_cases(path):
    """
    Loads extract_sequence test cases from a YAML file.

    Args:
        path (Path): The path to the YAML file containing the test cases.

    Returns:
        List[Any]: A list of pytest.param objects, each containing:
                   seq_id (str), win_start (int), win_end (int),
                   fasta_content (str), expected_sequence (str),
                   and a test_id.
                   Returns an empty list if the file is not found or is empty.
    """
    with open(path, "r", encoding="utf-8") as f:
        cases = yaml.safe_load(f)

    test_cases = []
    for case in cases:
        test_id = case.get("test_id")
        inputs = case.get("inputs", {})
        seq_id = inputs.get("seq_id")
        win_start = inputs.get("win_start")
        win_end = inputs.get("win_end")
        # Default fasta_content to an empty string if null in YAML
        fasta_content = (
            inputs.get("fasta_content")
            if inputs.get("fasta_content") is not None
            else ""
        )
        expected_sequence = case.get("expected_sequence")
        test_cases.append(
            pytest.param(
                seq_id, win_start, win_end, fasta_content, expected_sequence, id=test_id
            )
        )
    return test_cases
