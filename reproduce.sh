#!/usr/bin/env bash
set -euo pipefail

# Default gene list and chromosome mapping
declare -A gene_chr=(
    ["FANCA"]="16"
    ["FANCB"]="X"
    ["FANCC"]="9"
    ["BRCA2"]="13"
    ["FANCD2"]="3"
    ["FANCE"]="6"
    ["FANCF"]="11"
    ["FANCG"]="9"
    ["FANCI"]="15"
    ["BRIP1"]="17"
    ["FANCL"]="2"
    ["FANCM"]="14"
    ["PALB2"]="16"
    ["RAD51C"]="17"
    ["SLX4"]="16"
    ["ERCC4"]="16"
    ["RAD51"]="15"
    ["BRCA1"]="17"
    ["UBE2T"]="1"
    ["XRCC2"]="7"
    ["MAD2L2"]="1"
    ["RFWD3"]="16"
)

# Allow override of gene list via arguments
if [ "$#" -ge 1 ]; then
    GENES=("$@")
else
    GENES=("${!gene_chr[@]}")
fi

# Create directories
mkdir -p data data/scraped data/parsed data/BSData data/ClinVar_vcf data/LOVD_vcf

# Download and prepare ClinVar summary
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_GZ="data/variant_summary.txt.gz"
CLINVAR_TSV="data/variant_summary.tsv"
if [ ! -f "$CLINVAR_GZ" ]; then
    wget -O "$CLINVAR_GZ" "$CLINVAR_URL"
fi
if [ ! -f "$CLINVAR_TSV" ]; then
    gunzip -c "$CLINVAR_GZ" > "$CLINVAR_TSV"
fi

# Download and prepare reference genome
FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz"
FASTA_GZ="data/genome.fa.gz"
FASTA="data/genome.fa"
if [ ! -f "$FASTA" ]; then
    wget -O "$FASTA_GZ" "$FASTA_URL"
    gunzip -c "$FASTA_GZ" > "$FASTA"
    samtools faidx "$FASTA"
fi

# Process each gene
for gene in "${GENES[@]}"; do
    chr="${gene_chr[$gene]}"
    echo "Processing $gene on chromosome $chr"
    # ClinVar to VCF
    grep "GRCh37" "$CLINVAR_TSV" | grep -P "\\b${gene}\\b" | \
        cut -f 2,7,19,26,32,33,34 | \
        awk -F'\t' -v OFS='\t' '{print $3, $5, $6, $7, $1, $2, $4}' \
        > "data/ClinVar_vcf/${gene}.vcf"

    # LOVD scraping and parsing
    python3 get.py "$gene"
    python3 parse.py -L "data/scraped/${gene}_variants_full.tsv" --chr "$chr" -o "data/parsed/${gene}.tsv"
    python3 hgvsToVCF.py "GRCh37" "data/parsed/${gene}.tsv"
done
