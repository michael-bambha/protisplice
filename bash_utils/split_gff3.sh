#!/bin/bash

for chrom in $(cut -f1 gencode.v48.annotation.gff3 | grep -v '^#' | sort | uniq); do
  awk -v chr="$chrom" '$1 == chr || /^#/' gencode.v48.annotation.gff3 > "gencode.v48.annotation.${chrom}.gff3"
done
