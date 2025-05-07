#!/usr/bin/env python3
import argparse
import subprocess
import json
import sys

def main():
    parser = argparse.ArgumentParser(
        description="For each value in column 2 of a TSV, run mutalyzer_normalizer and print lines with errors")
    parser.add_argument('infile', help="Input TSV file")
    args = parser.parse_args()

    with open(args.infile, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 2:
                continue
            variant = cols[1]

            try:
                # invoke mutalyzer_normalizer command
                proc = subprocess.run(
                    ['mutalyzer_normalizer', variant],
                    capture_output=True, text=True, check=True
                )
                data = json.loads(proc.stdout)
            except Exception:
                # treat any failure as an error
                print(line)
                continue

            if 'ESYNTAXUC' in str(data):
                print(line)
            else:
                print(data)

if __name__ == '__main__':
    main()
