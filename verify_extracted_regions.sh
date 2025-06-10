#!/bin/bash

# Usage: ./verify_extracted_regions.sh annotations.gff3 out1.fasta

set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <annotations.gff3> <positive_seqs.fasta>"
    exit 1
fi

GFF_FILE="$1"
FASTA_FILE="$2"

EXONS_BED="exons.bed"
EXTRACTED_BED="extracted_windows.bed"


echo "Generating BED from extracted FASTA: $FASTA_FILE ..."
grep '^>' "$FASTA_FILE" | sed 's/^>//' | awk -F '_' 'NF >= 5 && $4 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ {
    OFS = "\t";
    print $1, $4 - 1, $5, $2, ".", $3;
}' > "$EXTRACTED_BED"

echo "Intersecting extracted regions with exons ..."
bedtools intersect -a extracted_windows.bed -b exons.bed -s -wo > matches_exon_overlap.bed

echo "Done. Output: matches_exon_overlap.bed"