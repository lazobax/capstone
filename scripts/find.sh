#!/bin/bash

gtf_file="$1"
coordinate="$2"
chrom_number="$3"

# Format chrom number to NC_0000XX
if [ "$chrom_number" -lt 10 ]; then
    formatted_chrom="NC_00000${chrom_number}"
else
    formatted_chrom="NC_0000${chrom_number}"
fi

awk -v coord="$coordinate" -v chrom="$formatted_chrom" '
BEGIN { FS=OFS="\t" }
{
    if ($1 ~ chrom && coord >= $4 && coord <= $5) {
        print $0
    }
}
' "$gtf_file"
