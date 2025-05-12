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



logfile="lines.log"
# Header: file, before, replaceCHROM, clean, bgzip, annotate, sort, index
echo -e "file\tbefore\treplaceCHROM\tclean\tbgzip\tannotate\tsort\tindex" > "$logfile"

for sub in LOVD CV; do
  for in_vcf in VCF_regular/${sub}/*.vcf; do
    base=$(basename "$in_vcf" .vcf)
    out_prefix="VCF/${sub}/${base}"

    # 1) count before
    cnt_before=$(grep -vc '^#' "$in_vcf")

    # a) replaceCHR → contig names
    python scripts/replaceCHROM.py "$in_vcf" "${out_prefix}.vcf"
    cnt_replace=$(grep -vc '^#' "${out_prefix}.vcf")

    # b) clean up IUPAC ambigs / non-ACGTN / negative POS / REF mismatches
    python scripts/removeIUPACambg.py \
      "${out_prefix}.vcf" \
      --clean "${out_prefix}.clean.vcf" \
      --error "${out_prefix}.error.vcf" \
      --mismatch "${out_prefix}.mismatch.vcf" \
      --ref data/genome.fa
    cnt_clean=$(grep -vc '^#' "${out_prefix}.clean.vcf")

    # c) bgzip the cleaned VCF
    bgzip -c "${out_prefix}.clean.vcf" > "${out_prefix}.vgz"
    cnt_bgzip=$(gunzip -c "${out_prefix}.vgz" | grep -vc '^#')

    # d) inject contig headers
    bcftools reheader \
      --fai data/genome.fa.fai \
      -o "${out_prefix}.contig.vcf.gz" \
      "${out_prefix}.vgz"
    cnt_annotate=$(gunzip -c "${out_prefix}.contig.vcf.gz" | grep -vc '^#')

    # e) sort & re-compress
    bcftools sort \
      "${out_prefix}.contig.vcf.gz" \
      -Oz -o "${out_prefix}.sorted.vcf.gz"
    cnt_sort=$(gunzip -c "${out_prefix}.sorted.vcf.gz" | grep -vc '^#')


    # f) index (final)
    tabix -p vcf "${out_prefix}.sorted.vcf.gz"
    cnt_index=$(bcftools view -H "${out_prefix}.sorted.vcf.gz" | wc -l)

    # append results to log
    echo -e "${sub}/${base}\t${cnt_before}\t${cnt_replace}\t${cnt_clean}\t${cnt_bgzip}\t${cnt_annotate}\t${cnt_sort}\t${cnt_index}" \
      >> "$logfile"
  done
done




# 4) Concatenate all sorted VCFs into one uncompressed VCF
bcftools concat -a -O v \
  VCF/LOVD/*.sorted.vcf.gz VCF/CV/*.sorted.vcf.gz \
  > data/combined.vcf 

echo "AFTER CONCAT"

echo $(cat data/combined.vcf | grep -cv "^#")


# removes non-ATCGN bases and keeps them, non ATCGN bases are IUPAC ambiguity
python scripts/removeIUPACambg.py



bgzip -c data/combined.clean.vcf  > data/combined.clean.vcf.gz
bcftools norm -m -both -f data/genome.fa \
  -Oz -cw -o data/combined.normalized.vcf.gz \
  data/combined.clean.vcf.gz
tabix -p vcf data/combined.normalized.vcf.gz

echo $(gunzip -c data/combined.normalized.vcf.gz | grep -cv "^#")



# Downloading ANNOVAR and the refGene, dbnsfp42a databases
# The following take some time. Especially dbnsfp42a, it's over 40Gb unzipped. 
# wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
# tar -xvzf annovar.latest.tar.gz
# perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene annovar/humandb/
# perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42a annovar/humandb/
# Rename both header and CHROM column without touching any records
gunzip -c data/combined.normalized.vcf.gz \
  | awk 'BEGIN{
      FS=OFS="\t";
      while(getline<"utils/acc2num.txt") map[$1]=$2
    }
    /^#/ { print; next }
    { if(map[$1]) $1=map[$1]; print }
  ' \
  | bgzip -c > data/combined.normalized.num.vcf.gz \
  && tabix -p vcf data/combined.normalized.num.vcf.gz

echo $(gunzip -c data/combined.normalized.num.vcf.gz | grep -cv "^#")

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

gunzip data/combined.normalized.num.vcf
mv data/combined.normalized.hg19_multianno.txt data/finalAnnot.tsv

mkdir -p data/combined/
mv data/combined.* data/combined/ 
mv data/combined/combined.normalized.vcf data/finalVCF.vcf

# 1) Download and prepare the chain file
mkdir -p data
wget -O data/19To38.chain.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip -f data/19To38.chain.gz

# 2) Download and index hg38 reference
wget -O data/hg38.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -f data/hg38.fa.gz
samtools faidx data/hg38.fa

# # 4) Liftover your VCF
CrossMap.py vcf \
  data/19To38.chain \
  data/finalVCF.num.vcf \
  data/hg38.fa \
  data/final.hg38.vcf

# #failed to map 14

bcftools reheader --fai data/hg38.fa.fai -o data/final.hg38.contig.vcf data/final.hg38.vcf

awk 'BEGIN{OFS="\t"} \
  /^#/ { print; next } \
  ($1 ~ /^[0-9]+$/ || $1=="X" || $1=="Y") { $1 = "chr"$1 } \
  { print }' \
  data/final.hg38.contig.vcf > data/final.hg38.nums.vcf

bgzip -c data/final.hg38.nums.vcf > data/final.hg38.nums.vcf.gz

bcftools norm -m -both -f data/hg38.fa -Oz -cw -o data/final.hg38.normalized.vcf.gz data/final.hg38.nums.vcf.gz 


