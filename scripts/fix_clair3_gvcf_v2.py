#!/usr/bin/env python

import argparse
import logging
import pathlib
import pysam
import sys
import gzip

def smart_open(fn, mode = "rt"):
    if str(fn).endswith(".gz"):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', action='store', dest='input_gvcf', type=pathlib.Path, required=True,
                        help='Input (zipped or unzipped) gVCF file path')
    parser.add_argument('-r', '--ref', action='store', dest='ref_fasta', type=pathlib.Path, required=True,
                        help='Reference FASTA file path')
    parser.add_argument('-o', '--output', action='store', dest='output_gvcf', type=pathlib.Path, required=True,
                        help='Output (unzipped) gVCF file path')

    return parser.parse_args()

def main():
    args = parse_arguments()

    if not pathlib.Path(f'{args.ref_fasta}.fai').exists():
        sys.exit(f'Index file does not exist for \'{args.ref_fasta}\'')

    ref_seq = pysam.FastaFile(args.ref_fasta)

    with smart_open(args.input_gvcf) as src, smart_open(args.output_gvcf, "w") as dst:
        for l in src:
            if l.startswith("chrY\t"):
                t = l.split('\t')
                if t[3] == "N":
                    pos = int(t[1])
                    t[3] = ref_seq.fetch("chrY", pos - 1, pos)
                    l = '\t'.join(t)
            dst.write(l)

    ref_seq.close()

if __name__ == '__main__':
    main()

