import sys
import pandas as pd
import hgvs
import hgvs.dataproviders.uta
from hgvs.exceptions import HGVSError
from bioutils.assemblies import make_ac_name_map, make_name_ac_map
from bioutils.sequences import reverse_complement
import hgvs.normalizer
from hgvs.edit import NARefAlt
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant
import os

def _as_interbase(posedit):
    if posedit.edit.type == "ins":
        start_i = posedit.pos.start.base
        end_i = posedit.pos.end.base - 1
    else:
        start_i = posedit.pos.start.base - 1
        end_i = posedit.pos.end.base
    return (start_i, end_i)

class Babelfish:
    def __init__(self, hdp, assembly_name):
        self.assembly_name = assembly_name
        self.hdp = hdp
        self.hn = hgvs.normalizer.Normalizer(
            hdp, cross_boundaries=False, shuffle_direction=5, validate=False
        )
        self.ac_to_name_map = make_ac_name_map(assembly_name)
        self.name_to_ac_map = make_name_ac_map(assembly_name)
        self.name_to_ac_map.update({ac: ac for ac in self.name_to_ac_map.values()})

    def hgvs_to_vcf(self, var_g):
        if var_g.type != "g":
            raise RuntimeError(f"Expected g. variant, got {var_g}")
        vleft = self.hn.normalize(var_g)
        (start_i, end_i) = _as_interbase(vleft.posedit)
        chrom = self.ac_to_name_map[vleft.ac]
        typ = vleft.posedit.edit.type

        if typ == "dup":
            start_i -= 1
            alt = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
            ref = alt[0]
        elif typ == "inv":
            ref = vleft.posedit.edit.ref
            alt = reverse_complement(ref)
        else:
            alt = vleft.posedit.edit.alt or ""
            if typ in ("del", "ins"):
                start_i -= 1
                ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
                alt = ref[0] + alt
            else:
                ref = vleft.posedit.edit.ref
                if ref == alt:
                    alt = "."
        return chrom, start_i + 1, ref, alt, typ

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <assembly_name> <input.tsv>")
        sys.exit(1)
    assembly_name = sys.argv[1]
    input_tsv = sys.argv[2]
    output_vcf = "data/LOVD_vcf/" + os.path.basename(input_tsv).rsplit(".", 1)[0] + ".vcf"


    hdp = hgvs.dataproviders.uta.connect()
    import pysam
    from functools import lru_cache

    # 1) open your local FASTA once
    fasta = pysam.FastaFile("data/genome.fa")

    # 2) cache every fetch to avoid repeated I/O
    @lru_cache(maxsize=None)
    def fetch_local(ac, start, end):
        # pysam.fetch uses 0-based, half-open intervals
        return fasta.fetch(ac, start, end)

    # 3) monkey-patch the UTA seqfetcher
    hdp = hgvs.dataproviders.uta.connect()
    hdp.seqfetcher.fetch_seq = fetch_local

    bf = Babelfish(hdp, assembly_name)
    hp = hgvs.parser.Parser()


    df = pd.read_csv(input_tsv, sep="\t", header=None)
    failed = []
    vcf_entries = []

    for idx, row in df.iterrows():
        hgvs_str = str(row[1]).strip()
        ID = str(row[0]).strip()
        CONSEQ = str(row[2]).strip()
        REPORTED = str(row[3]).strip()
        try:
            var_g = hp.parse_hgvs_variant(hgvs_str)
            chrom, pos, ref, alt, typ = bf.hgvs_to_vcf(var_g)
            # print(chrom, len(pos), ref, alt, typ )
            vcf_entries.append([chrom, pos, ref, alt,typ, ID, CONSEQ, REPORTED])
        except Exception as e:
            failed.append(hgvs_str)

    vcf_df = pd.DataFrame(vcf_entries, columns=["#CHROM", "POS", "REF", "ALT", "TYP", "ID", "CONSEQ", "REPORTED"])
    vcf_df.to_csv(output_vcf, sep="\t", index=False)

    if failed:
        print("Failed to convert the following variants:")
        for f in failed:
            print(f)
        # with open(f"{os.path.basename(input_tsv).rsplit(".", 1)[0]}_failedConversion.txt", "w") as file:
        #     file.write("\n".join(failed))
    else:
        print("All variants converted successfully!")
