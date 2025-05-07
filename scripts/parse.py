#!/usr/bin/env python3
import argparse
import csv
import re
import sys
import os

def format_contig(chrom):
    chrDict = { "1":"NC_000001.10",
                "2":"NC_000002.11",
                "3":"NC_000003.11",
                "4":"NC_000004.11",
                "5":"NC_000005.9",
                "6":"NC_000006.11",
                "7":"NC_000007.13",
                "8":"NC_000008.10",
                "9":"NC_000009.11",
                "10":"NC_000010.10",
                "11":"NC_000011.9",
                "12":"NC_000012.11",
                "13":"NC_000013.10",
                "14":"NC_000014.8",
                "15":"NC_000015.9",
                "16":"NC_000016.9",
                "17":"NC_000017.10",
                "18":"NC_000018.9",
                "19":"NC_000019.9",
                "20":"NC_000020.10",
                "21":"NC_000021.8",
                "22":"NC_000022.10",
                "X":"NC_000023.10",
                "Y":"NC_000024.9"}
    return chrDict.get(chrom, "")


def process_c_file(path, contig, ignored_writer):
    """
    Process the -C TSV:
      - skip header
      - skip rows where col2 split on '|' has >2 parts (log to ignored)
      - extract transcript + c.change from col1, else log to ignored
      - assign CV##### IDs
    Returns a list of (ID, VAR).
    """
    results = []
    counter = 1
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)  # skip header
        for row in reader:
            # filter on column 2
            if len(row) < 2 or len(row[1].split('|')) > 3:
                ignored_writer.writerow(row)
                continue

            # parse column 1 for transcript and c.change
            m = re.match(r'^([^(:]+)[^:]*:(c\.[^ )]+)', row[0])
            if not m:
                ignored_writer.writerow(row)
                continue
            transcript, change = m.group(1), m.group(2)

            var = f"{contig}({transcript}):{change}"
            vid = f"CV{counter:05d}"
            results.append((vid, var))
            counter += 1
    return results


def process_l_file(path, contig, ignored_writer, gene):
    """
    Process the -L TSV:
      - skip header
      - take the 'g.…' string from col11, else log to ignored
      - prefix with the given contig
      - assign LOVD##### IDs
    """
    results = []
    counter = 1
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 11:
                ignored_writer.writerow(row)
                continue
            gpart = row[9].strip()
            if not gpart.startswith('g.'):
                ignored_writer.writerow(row)
                continue
            conseq = row[8].strip()
            reported = row[1].strip()
            var = f"{contig}:{gpart}"
            vid = f"LOVD.{gene}{counter:05d}"
            results.append((vid, var, conseq, reported))
            counter += 1
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Merge -C and -L TSVs into ID<TAB>VAR")
    parser.add_argument('-C', '--cfile', required=False, help="input TSV for -C")
    parser.add_argument('-L', '--lfile', required=True, help="input TSV for -L")
    parser.add_argument('-o', '--outfile', required=True, help="output TSV file")
    parser.add_argument('--chr', required=True, help="chromosome number (e.g. 16)")
    args = parser.parse_args()

    contig = format_contig(args.chr)
    gene = os.path.basename(args.lfile).split('_')[0]
    # open ignored rows file
    ignored_path = 'data/ignored/' + gene + '_ignored.tsv'
    with open(ignored_path, 'w', newline='') as ign_f, \
         open(args.outfile, 'w', newline='') as out_f:
        ignored_writer = csv.writer(ign_f, delimiter='\t', lineterminator='\n')
        out_writer = csv.writer(out_f, delimiter='\t', lineterminator='\n')

        # write header for main output
        out_writer.writerow(['ID', 'VAR', 'CONSEQ', 'REPORTED'])

        # process C and write
        # c_entries = process_c_file(args.cfile, contig, ignored_writer)
        # for vid, var in c_entries:
        #     out_writer.writerow([vid, var])

        # process L and write
        l_entries = process_l_file(args.lfile, contig, ignored_writer, gene)
        for vid, var, conseq, reported in l_entries:
            out_writer.writerow([vid, var, conseq, reported])

    print(f"Merged variants written to {args.outfile}")
    print(f"Ignored rows written to {ignored_path}")

if __name__ == '__main__':
    main()



