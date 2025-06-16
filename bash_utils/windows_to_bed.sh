#!/bin/bash

grep '^>' | sed 's/>//' | awk -F '_' '{
    print $1"\t"($5-1)"\t"$6"\t"$2"\t.\t"$3;
}' > extracted_windows.bed

