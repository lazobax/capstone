#!/usr/bin/env python3
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Replace chromosome names in a VCF using a mapping dictionary"
    )
    parser.add_argument(
        "input_vcf",
        help="Path to the input VCF file"
    )
    parser.add_argument(
        "output_vcf",
        help="Path to write the updated VCF file"
    )
    return parser.parse_args()

# Mapping from simple chromosome labels to contig accessions
CHR_DICT = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "X": "NC_000023.10",
    "Y": "NC_000024.9"
}


def main():
    args = parse_args()

    try:
        infile = open(args.input_vcf, 'r')
    except Exception as e:
        sys.exit(f"Error opening input VCF: {e}")

    try:
        outfile = open(args.output_vcf, 'w')
    except Exception as e:
        sys.exit(f"Error opening output VCF: {e}")

    for line in infile:
        if line.startswith('#'):
            # Preserve header lines unchanged
            outfile.write(line)
            continue

        fields = line.rstrip('\n').split('\t')
        chrom = fields[0]
        # Replace if in mapping
        if chrom in CHR_DICT:
            fields[0] = CHR_DICT[chrom]
        # else leave as-is
        outfile.write("\t".join(fields) + "\n")

    infile.close()
    outfile.close()

if __name__ == '__main__':
    main()
