#!/bin/bash

grep '\texon\t' annotations.gff3 | awk -F'\t' '{
    split($9, attrs, ";");
    for (i in attrs) {
        if (attrs[i] ~ /^Parent=/) {
            split(attrs[i], id, "=");
            transcript_id = id[2];
        }
    }
    print $1"\t"($4-1)"\t"$5"\t"transcript_id"\t.\t"$7;
}' > exons.bed