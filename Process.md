# Pipeline 

declare the genes and their respective chromosomes

```
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

mkdir data
mkdir data/scraped data/ClinVar_vcf/
mkdir data/parsed/
mkdir data/LOVD_vcf

```

## Get ClinVar datasets. 
for each gene do:
```
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz
grep "GRCh37" data/variant_summary.tsv |\
                        grep "${gene}\b" |\
             cut -f 2,7,19,26,32,33,34 |\
             awk -F'\t' -v OFS='\t' '{print $3, $5, $6, $7, $1, $2, $4}'\
             > data/ClinVar_vcf/${gene}.vcf
```

## Scrape LOVD data
for each gene do:

```
python3 get.py ${gene}
python scripts/parse.py -L data/scraped/${gene}_variants_full.tsv --chr ${chr_num} -o data/parsed/${gene}.tsv
python scripts/hgvsToVCF.py "GRCh37" data/parsed/${gene}.tsv  
```