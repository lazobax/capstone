#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <reference.fa> <input.vcf> <output.vcf>"
  exit 1
fi

ref="$1"
in_vcf="$2"
out_vcf="$3"
threads=6

# 1) Index reference if needed
if [ ! -f "${ref}.fai" ]; then
  echo "[+] Indexing reference: ${ref}"
  samtools faidx "$ref"
fi

# 2) Build a meta‐header with existing ##… plus our contigs
bcftools view -h "$in_vcf" | grep '^##' > header.meta
cut -f1,2 "${ref}.fai" \
  | awk '{ printf("##contig=<ID=%s,length=%s>\n",$1,$2) }' \
  >> header.meta
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> header.meta

# 3) Re‐header
echo "[+] Adding contigs to VCF header"
bcftools reheader -h header.meta "$in_vcf" > withcontigs.vcf

# 4) Filter out illegal POS (<=0)
echo "[+] Dropping any POS < 1"
awk '$1~/^#/ || $2>0' withcontigs.vcf > clean.vcf

# extract only biallelic SNPs
bcftools view \
  -v snps \
  clean.vcf \
  -Ov -o only_snps.vcf

# 5) Left‐normalize indels
echo "[+] Normalizing (left‐aligning) against ${ref}"
bcftools norm \
  -f "$ref" \
  --threads "$threads" \
  -c x \
  -Ov only_snps.vcf \
  > "$out_vcf"

# 6) Cleanup
rm -f header.meta withcontigs.vcf clean.vcf

echo "[✓] Done – output is $out_vcf"
