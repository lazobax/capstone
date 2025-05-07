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
mkdir -p data data/scraped data/parsed data/BSData data/ClinVar_vcf data/LOVD_vcf data/ignored/

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
    python3 scripts/addClinVarID.py ${gene}.vcf

    # LOVD scraping and parsing
    python3 scripts/get.py "$gene"
    python3 scripts/fixScrape.py "$gene"_variants_full.tsv
    python3 scripts/parse.py -L "data/scraped_adj/${gene}_variants_full.tsv" --chr "$chr" -o "data/parsed/${gene}.tsv"
    python3 scripts/hgvsToVCF.py "GRCh37" "data/parsed/${gene}.tsv"
done

# Reformating into valid, actual, VCF format
mkdir -p VCF_regular/CV VCF_regular/LOVD
python3 scripts/tsvToVcf.py --indir data/ClinVar_vcf/ --outdir VCF_regular/CV/
python3 scripts/tsvToVcf.py --indir data/LOVD_vcf/ --outdir VCF_regular/LOVD/


# 1) Build contigs.txt from genome.fa.fai
awk '{ printf("##contig=<ID=%s,length=%s>\n",$1,$2) }' \
    data/genome.fa.fai > data/contigs.txt

# 2) Prepare output directories
mkdir -p VCF/LOVD VCF/CV

# 3) Process every VCF in VCF_regular/LOVD and VCF_regular/CV
for sub in LOVD CV; do
  for in_vcf in VCF_regular/${sub}/*.vcf; do
    base=$(basename "$in_vcf" .vcf)
    out_vcf=VCF/${sub}/${base}.vcf

    # a) Replace chr → contig names
    python scripts/replaceCHROM.py "$in_vcf" "$out_vcf"

    # b) Compress with bgzip
    bgzip -c "$out_vcf" > "${out_vcf}.gz"

    # c) Inject contig headers
    bcftools annotate \
      --header-lines data/contigs.txt \
      --rename-chrs utils/num2acc.txt \
      -Oz -o "${out_vcf%.vcf}.contig.vcf.gz" \
      "${out_vcf}.gz"

    # d) Sort and re-compress
    bcftools sort \
      "${out_vcf%.vcf}.contig.vcf.gz" \
      -Oz -o "${out_vcf%.vcf}.sorted.vcf.gz"

    # e) Index
    tabix -p vcf "${out_vcf%.vcf}.sorted.vcf.gz"
  done
done

# 4) Concatenate all sorted VCFs into one uncompressed VCF
bcftools concat -a -O v \
  VCF/LOVD/*.sorted.vcf.gz VCF/CV/*.sorted.vcf.gz \
  > data/combined.vcf

#removes non-ATCGN bases and keeps them, non ATCGN bases are IUPAC ambiguity
python scripts/removeIUPACambg.py

# 1) bgzip & index the raw combined VCF
bgzip -c data/combined.clean.vcf   > data/combined.clean.vcf.gz
tabix -p vcf data/combined.clean.vcf.gz
echo "zipped"

# 2) Split multi-allelic sites into biallelic records
bcftools norm -m -both \
  data/combined.clean.vcf.gz \
  -Oz -o data/combined.split.vcf.gz
echo "normed"

# 3) Left-align and trim indels against your reference
bcftools norm -f data/genome.fa \
  data/combined.split.vcf.gz \
  -Oz -cw -o data/combined.normalized.vcf.gz
echo "normed 2" 

# 4) Index the normalized VCF
tabix -p vcf data/combined.normalized.vcf.gz



# Downloading ANNOVAR and the refGene, dbnsfp42a databases
# The following take some time. Especially dbnsfp42a, it's over 40Gb unzipped. 
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -xvzf annovar.latest.tar.gz
perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene annovar/humandb/
perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42a annovar/humandb/

bcftools annotate \
      --rename-chrs utils/acc2num.txt \
      -Oz -o data/combined.normalized.num.vcf.gz \
      data/combined.normalized.vcf.gz


perl annovar/table_annovar.pl \
    data/combined.normalized.num.vcf.gz \
    annovar/humandb/ \
    -buildver hg19 \
    -out data/combined.normalized \
    -remove \
    -protocol refGene,dbnsfp42a \
    -operation g,f \
    -nastring . \
    -vcfinput

mv data/combined.normalized.hg19_multianno.txt data/finalAnnot.tsv

mkdir -p data/combined/
mv data/combined.* data/combined/ 