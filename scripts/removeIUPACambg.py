#!/usr/bin/env python3
import os
import pysam

def main():
    input_vcf = 'data/combined.vcf'
    clean_vcf = 'data/combined.clean.vcf'
    error_vcf = 'data/error.vcf'
    mismatch_vcf = 'data/ref_mismatch.vcf'

    # Ensure output directories exist
    for path in (clean_vcf, error_vcf, mismatch_vcf):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # Allowed bases for REF/ALT
    allowed = set('ACGTN')
    # Load reference genome
    fasta = pysam.FastaFile('data/genome.fa')

    with open(input_vcf, 'r') as fin, \
         open(clean_vcf, 'w') as fout_good, \
         open(error_vcf, 'w') as fout_err, \
         open(mismatch_vcf, 'w') as fout_mis:

        for line in fin:
            if line.startswith('#'):
                # Preserve headers in all outputs
                fout_good.write(line)
                fout_err.write(line)
                fout_mis.write(line)
                continue

            cols = line.rstrip('\n').split('\t')
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alts = cols[4].split(',')

            # Check for standard bases
            if not (all(b in allowed for b in ref) and all(all(b in allowed for b in alt) for alt in alts)):
                fout_err.write(line)
                continue

            # Validate REF matches reference genome
            start = pos - 1
            end = start + len(ref)
            try:
                seq = fasta.fetch(chrom, start, end)
            except Exception:
                fout_mis.write(line)
                continue

            if seq.upper() != ref.upper():
                fout_mis.write(line)
            else:
                fout_good.write(line)

    print(f"Clean VCF written to: {clean_vcf}")
    print(f"Non-standard records in: {error_vcf}")
    print(f"Reference mismatches in: {mismatch_vcf}")

if __name__ == '__main__':
    main()
