#!/usr/bin/env python

import argparse
import logging
import pathlib
import pysam
import sys
import vcf

def create_log_file_path(input_gvcf, output_gvcf):
    return str(output_gvcf.with_name(f'{input_gvcf.name}.ref_base_fix.log'))

def append_missing_non_ref_allele_freq(vcf_record, sample_id):
    sample_var_call = vcf_record.genotype(sample_id)
    sample_var_call.data.AF.append(sample_var_call.data.AD[-1] / sample_var_call.data.DP)

    return vcf_record

def correct_ref_base(vcf_record, ref_seq, chrom='chrY'):
    ref_base = ref_seq.fetch(chrom, vcf_record.POS - 1, vcf_record.POS)

    if ref_base != str(vcf_record.REF):
        logging.info(f'{chrom}:{vcf_record.POS} - Change REF from {str(vcf_record.REF)} to {ref_base}')
        vcf_record.REF = ref_base

    return vcf_record

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

    gvcf_reader = vcf.Reader(filename=str(args.input_gvcf))

    if len(gvcf_reader.samples) > 1:
        sys.exit(f'The gVCF file \'{args.input_gvcf}\' must contain one sample only')

    sample_id = gvcf_reader.samples[0]

    logging.basicConfig(
        filename=create_log_file_path(args.input_gvcf, args.output_gvcf), level=logging.INFO, format='%(message)s'
    )

    with open(args.output_gvcf, 'w') as gvcf_out:
        gvcf_writer = vcf.Writer(gvcf_out, gvcf_reader)

        for vcf_record in gvcf_reader:
            if str(vcf_record.ALT[-1]) == '<NON_REF>':
                if len(vcf_record.ALT) > 1:
                    gvcf_writer.write_record(append_missing_non_ref_allele_freq(vcf_record, sample_id))
                    continue
                elif str(vcf_record.CHROM) == 'chrY' and str(vcf_record.REF) == 'N':
                    gvcf_writer.write_record(correct_ref_base(vcf_record, ref_seq))
                    continue
            
            gvcf_writer.write_record(vcf_record)

        gvcf_writer.close()

    ref_seq.close()

if __name__ == '__main__':
    main()
