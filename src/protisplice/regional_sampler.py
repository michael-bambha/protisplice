"""
File: regional_sampler.py
Author: Michael Bambha
Description: Module for sampling genomic regions (intergenic,
exonic, and intronic).
"""

# pylint:disable=no-member

import random
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from collections import defaultdict
import pysam

from .data_models import (
    SpliceJunction,
    JunctionData,
    StrandType,
    JunctionType,
    Transcript,
    Gene,
    SamplingParams,
)

from .sequence_extraction import extract_sequence, add_intron_coords
from .utils import apply_strand


class RegionalSampler:
    """Randomly sample exons, introns, and intergenic regions"""

    def __init__(self, fasta_path: str, sampling_params: SamplingParams):
        self.fasta_path = Path(fasta_path)
        self.params = sampling_params

        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}.")

        fai_path = Path(f"{self.fasta_path}.fai")

        if not fai_path.exists():
            print(f"FASTA index {fasta_path} not found -- generating index...")
            pysam.faidx(fasta_path)
            print(f"FASTA index {fasta_path}.fai generated.")

    def sample_introns(
        self, transcripts: Dict[str, Transcript], target_count: int, seed: int = 100
    ) -> List[JunctionData]:
        """Sample introns from transcripts.
        Args:
            transcripts (Dict[str, Transcript]): Dict of transcript_id: Transcript
            target_count (int): Max number of introns to sample
            seed (int, optional): Random state. Defaults to 100.

        Returns:
            List[JunctionData]: List of JunctionData objects in the format
            {junction: SpliceJunction, win_start: win_start, win_end: win_end, seq: seq}.
        """
        random.seed(seed)
        sequences = []

        # Add intron coordinates to transcripts
        transcripts_with_introns = add_intron_coords(transcripts)

        with pysam.FastaFile(str(self.fasta_path)) as fasta:
            transcript_ids = list(transcripts_with_introns.keys())
            random.shuffle(transcript_ids)

            for transcript_id in transcript_ids:
                if len(sequences) >= target_count:
                    break

                transcript = transcripts_with_introns[transcript_id]
                extracted = self._sample_from_transcript(
                    fasta, transcript_id, transcript, target_count - len(sequences)
                )
                sequences.extend(extracted)

        return sequences

    def sample_exonic_regions(
        self, transcripts: Dict[str, Transcript], target_count: int, seed: int = 100
    ) -> List[JunctionData]:
        """Sample exons from a dictionary of transcripts.

        Args:
            transcripts (Dict[str, Transcript]): Dict of transcript_id: Transcript
            target_count (int): Max number of exons to sample
            seed (int, optional): Random state. Defaults to 100.

        Returns:
            List[JunctionData]: List of JunctionData objects in the format
            {junction: SpliceJunction, win_start: win_start, win_end: win_end, seq: seq}.
        """
        random.seed(seed)
        sequences = []

        exons = []
        for transcript_id, transcript in transcripts.items():
            for i, (exon_start, exon_end) in enumerate(transcript.exons):
                exons.append(
                    {
                        "transcript_id": transcript_id,
                        "exon_index": i,
                        "start": exon_start,
                        "end": exon_end,
                        "seqid": transcript.info.seqid,
                        "strand": transcript.info.strand,
                    }
                )

        random.shuffle(exons)

        with pysam.FastaFile(str(self.fasta_path)) as fasta:
            for exon_data in exons:
                if len(sequences) >= target_count:
                    break
                junction_data = self._sample_from_exon(fasta, exon_data, len(sequences))
                if junction_data:
                    sequences.append(junction_data)

        return sequences

    def sample_intergenic_regions(
        self,
        genes: Dict[str, Gene],
        chromosome_lengths: Dict[str, int],
        target_count: int,
        seed: int = 100,
    ) -> List[JunctionData]:
        """Sample regions between genes, given a dictionary of genes,
        and length of the chromosomes.

        Args:
            genes (Dict[str, Gene]): Dict of gene_id: Gene
            chromosome_lengths (Dict[str, int]): Dict of chromosome ID: length
            target_count (int): Number of sequences to sample
            seed (int, optional): Random state. Defaults to 100.

        Returns:
            List[JunctionData]: List of JunctionData objects in the format
            {junction: SpliceJunction, win_start: win_start, win_end: win_end, seq: seq}.
        """
        random.seed(seed)
        sequences = []

        intergenic_regions = self._identify_intergenic_regions(
            genes, chromosome_lengths
        )

        with pysam.FastaFile(str(self.fasta_path)) as fasta:
            region_list = []
            for seq_id, regions in intergenic_regions.items():
                region_list.extend([(seq_id, region) for region in regions])

            random.shuffle(region_list)

            for seq_id, region in region_list:
                if len(sequences) >= target_count:
                    break

                junction_data = self._sample_from_intergenic_region(
                    fasta, seq_id, region, len(sequences)
                )
                if junction_data:
                    sequences.append(junction_data)

        return sequences

    def _sample_from_transcript(
        self,
        fasta: pysam.FastaFile,
        transcript_id: str,
        transcript: Transcript,
        max_samples: int,
    ) -> List[JunctionData]:
        """Sample introns of a transcript. Can also provide a maximum cap on the number of samples,
        but if downstream buffer_size is moderate to high, then max_samples will not be reached.

        Args:
            fasta (pysam.FastaFile): pysam FastaFile object
            transcript_id (str): ID of transcript
            transcript (Transcript): Transcript data class
            max_samples (int): maximum cap on # samples to obtain

        Returns:
            List[JunctionData]: List of JunctionData objects in the format
            {junction: SpliceJunction, win_start: win_start, win_end: win_end, seq: seq}.
        """
        if not transcript.introns:
            return []

        sequences = []

        for intron_start, intron_end in transcript.introns:
            if len(sequences) >= max_samples:
                break

            junction_data = self._sample_from_intron(
                fasta,
                transcript_id,
                transcript,
                intron_start,
                intron_end,
                len(sequences),
            )
            if junction_data:
                sequences.append(junction_data)

        return sequences

    def _sample_from_intron(
        self,
        fasta: pysam.FastaFile,
        transcript_id: str,
        transcript: Transcript,
        intron_start: int,
        intron_end: int,
        sample_index: int,
    ) -> Optional[JunctionData]:
        """Sample sequence regions from a singular intron with defined start/end coordinates.

        Args:
            fasta (pysam.FastaFile): pysam FASTA object
            transcript_id (str): ID of transcript
            transcript (Transcript): Transcript data class
            intron_start (int): Desired start coordinate
            intron_end (int): Desired end coordinate
            sample_index (int): Number of the intron ordered 5' to 3'.

        Returns:
            Optional[JunctionData]: List of JunctionData objects in the format
            {junction: SpliceJunction, win_start: win_start, win_end: win_end, seq: seq}.
        """
        rand_start, rand_end = self._random_window(
            intron_start, intron_end, self.params.window_size, self.params.buffer_size
        )
        if rand_start is None or rand_end is None:
            return None

        seq = extract_sequence(fasta, transcript.info.seqid, rand_start, rand_end)
        if not seq:
            return None

        strand = transcript.info.strand
        seq = apply_strand(seq, strand)

        pseudo_junction = self._create_pseudo_junction(
            seqid=transcript_id,
            coord=rand_start,
            strand=transcript.info.strand,
            junction_type=JunctionType.INTRON,
            sample_index=sample_index + 1,
        )

        return JunctionData(
            junction=pseudo_junction,
            window_start=rand_start,
            window_end=rand_end,
            sequence=seq,
        )

    def _sample_from_exon(
        self, fasta: pysam.FastaFile, exon_data: Dict, sample_index: int
    ) -> Optional[JunctionData]:
        """Sample a sequence from a random start position in an exon.

        Args:
            fasta (pysam.FastaFile): pysam FASTA object
            exon_data (dict): {
                        "transcript_id": transcript_id,
                        "exon_index": i,
                        "start": exon_start,
                        "end": exon_end,
                        "seqid": transcript.info.seqid,
                        "strand": transcript.info.strand,
                    }
            sample_index (int): Number of the exon from 5' to 3'

        Returns:
            Optional[JunctionData]: JunctionData object in the format
            {junction: SpliceJunction, win_start: win_start, win_end: win_end, seq: seq}.
        """
        exon_start = exon_data["start"]
        exon_end = exon_data["end"]

        rand_start, rand_end = self._random_window(
            exon_start, exon_end, self.params.window_size, self.params.buffer_size
        )
        if rand_start is None or rand_end is None:
            return None

        seq = extract_sequence(fasta, exon_data["seqid"], rand_start, rand_end)
        if not seq:
            return None

        strand = exon_data["strand"]
        seq = apply_strand(seq, strand)

        pseudo_junction = self._create_pseudo_junction(
            seqid=exon_data["transcript_id"],
            coord=rand_start,
            strand=strand,
            junction_type=JunctionType.EXON,
            sample_index=sample_index,
        )

        return JunctionData(
            junction=pseudo_junction,
            window_start=rand_start,
            window_end=rand_end,
            sequence=seq,
        )

    def _sample_from_intergenic_region(
        self,
        fasta: pysam.FastaFile,
        seqid: str,
        region: Tuple[int, int],
        sample_index: int,
    ) -> Optional[JunctionData]:
        """Given known coordinates, sample a sequence
        from a random start position in an intergenic region.

        Args:
            fasta (pysam.FastaFile): pysam FASTA object
            seqid (str): Name of the chromosome from FASTA header
            region (Tuple[int, int]): Tuple of (start, end) coords
            sample_index (int): Number of the intergenic region

        Returns:
            Optional[JunctionData]: JunctionData object in the format
            {junction: SpliceJunction, win_start: win_start, win_end: win_end, seq: seq}.
        """
        region_start, region_end = region

        rand_start, rand_end = self._random_window(
            region_start, region_end, self.params.window_size, self.params.buffer_size
        )
        if rand_start is None or rand_end is None:
            return None

        seq = extract_sequence(fasta, seqid, rand_start, rand_end)
        if not seq:
            return None

        # For intergenic regions, randomly choose strand
        strand = random.choice([StrandType.POSITIVE, StrandType.NEGATIVE])
        seq = apply_strand(seq, strand)

        pseudo_junction = self._create_pseudo_junction(
            seqid=seqid,
            coord=rand_start,
            strand=strand,
            junction_type=JunctionType.INTERGENIC,
            sample_index=sample_index,
        )

        return JunctionData(
            junction=pseudo_junction,
            window_start=rand_start,
            window_end=rand_end,
            sequence=seq,
        )

    def _identify_intergenic_regions(
        self,
        genes: Dict[str, Gene],
        chromosome_lengths: Dict[str, int],
    ) -> Dict[str, List[Tuple[int, int]]]:
        """Find intergenic regions, given a dictionary of genes
        and the length of each chromosome.

        Args:
            genes (Dict[str, Gene]): Dict of gene ID: Gene object
            chromosome_lengths (Dict[str, int]): Dict of chromosome #: length

        Returns:
            Dict[str, List[Tuple[int, int]]]: Dict of chromosome #: List[Tuple(
            start, end)] where start, end are the start and end coordinates
            of the intergenic region.
        """
        intergenic_regions = {}

        # group genes by chromosome
        genes_by_chr = defaultdict(list)
        for gene in genes.values():
            genes_by_chr[gene.seq_id].append(gene)

        for seqid, chr_genes in genes_by_chr.items():
            if seqid not in chromosome_lengths:
                continue

            sorted_genes = sorted(chr_genes, key=lambda x: x.start)
            regions = []

            # add region before first gene
            if sorted_genes and sorted_genes[0].start > 1:
                regions.append((1, sorted_genes[0].start - 1))

            # add regions between genes
            for i in range(len(sorted_genes) - 1):
                current_gene = sorted_genes[i]
                next_gene = sorted_genes[i + 1]

                if next_gene.start > current_gene.end + 1:
                    regions.append((current_gene.end + 1, next_gene.start - 1))

            # add region after last gene
            if sorted_genes and sorted_genes[-1].end < chromosome_lengths[seqid]:
                regions.append((sorted_genes[-1].end + 1, chromosome_lengths[seqid]))

            intergenic_regions[seqid] = regions

        return intergenic_regions

    def _random_window(
        self, region_start: int, region_end: int, window_size: int, buffer_size: int
    ) -> Tuple[Optional[int], Optional[int]]:
        """Given the window size, start and end coordinates of a region, and a buffer size, randomly
        choose valid window coordinates.

        Args:
            region_start (int): Coordinate of gene/exon/intron start
            region_end (int): Coordinate of gene/exon/intron end
            window_size (int): Number of bases for the sequence window
            buffer_size (int): Buffer size away from region_start/region_end

        Returns:
            Tuple[int, int]: Coordinates of (start, end) of the random window
        """
        safe_start = region_start + buffer_size
        safe_end = region_end - buffer_size

        if safe_end - safe_start + 1 < window_size:
            return None, None

        max_start = safe_end - window_size + 1
        if safe_start > max_start:
            return None, None

        rand_start = random.randint(safe_start, max_start)
        rand_end = rand_start + window_size - 1

        return rand_start, rand_end

    def _create_pseudo_junction(
        self,
        seqid: str,
        coord: int,
        strand: StrandType,
        junction_type: JunctionType,
        sample_index: int,
    ) -> SpliceJunction:
        """Create pseudo-junction for decoy sites to keep data format the same across decoy
        and true splice sites. Coord for pseudo-junctions should just be the start position
        of the random window (rand_start).

        Args:
            seqid (str): ID of the sampled junction
            coord (int): rand_start of the window
            strand (StrandType): (+) or (-)
            junction_type (JunctionType): "intergenic", "exon", or "intron"

        Returns:
            SpliceJunction: SpliceJunction object with
        """
        return SpliceJunction(
            id=f"{seqid}_{junction_type.value}_{sample_index}",
            seqid=seqid,
            coord=coord,
            strand=strand,
            junction_type=junction_type,
        )

    def get_chromosome_lengths(self) -> Dict[str, int]:
        """Find the length of any FASTA sequence with pysam.
        Implemented for finding chromosome lengths.

        Returns:
            Dict[str, int]: Dict of seq_id: length
        """
        lengths = {}
        with pysam.FastaFile(str(self.fasta_path)) as fasta:
            for ref in fasta.references:
                lengths[ref] = fasta.get_reference_length(ref)
        return lengths
