import os
import glob

vcf_header = """##fileformat=VCFv4.2
##INFO=<ID=ID,Number=1,Type=String,Description="Variant ID">
##INFO=<ID=CONSEQ,Number=1,Type=String,Description="Clinical consequence">
##INFO=<ID=REPORTED,Number=1,Type=Integer,Description="Number of times variant reported">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""

def convert_tsv_to_vcf(tsv_path, out_dir):
    basename = os.path.basename(tsv_path).replace(".vcf", "")
    vcf_path = os.path.join(out_dir, basename + ".vcf")

    with open(tsv_path) as tsv, open(vcf_path, "w") as vcf:
        vcf.write(vcf_header)
        next(tsv)  # skip header
        for line in tsv:
            cols = line.strip().split("\t")
            chrom, pos, ref, alt, typ, var_id, conseq, reported = cols
            info = f"ID={var_id};CONSEQ={conseq};REPORTED={reported}"
            vcf.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\t.\t{info}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", required=True, help="Input directory of TSV files")
    parser.add_argument("--outdir", required=True, help="Output directory for VCF files")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    for tsv_file in glob.glob(os.path.join(args.indir, "*.vcf")):
        convert_tsv_to_vcf(tsv_file, args.outdir)
