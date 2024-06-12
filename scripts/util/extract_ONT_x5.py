#!/usr/bin/env python3

from argparse import ArgumentParser
from os import cpu_count
from pathlib import Path, PurePosixPath
import numpy as np
import pod5
import pysam
import random
import re
import shlex
import subprocess
import sys

ONT_FAST5_SUBSET_TOOL_PATH = '/group/bioi1/jimmyc/conda/envs/python39/lib/python3.9/site-packages/ont_fast5_api/conversion_tools/fast5_subset.py'

FAST5_FILE_EXT = '.fast5'
POD5_FILE_EXT = '.pod5'

def is_input_chrom_valid(chrom):
    return re.match(r'^(chr)?([1-9]|1[0-9]{1}|2[012]{1}|[XYM])$', chrom) is not None

def is_user_label_valid(user_label):
    return re.match(r'^[\w\-\.]+$', user_label) is not None

def is_subsample_val_float(subsample_val):
    return re.match(r'^0\.[0-9]+$', subsample_val) is not None

def is_subsample_val_int(subsample_val):
    return re.match(r'^[0-9]+$', subsample_val) is not None

def is_subsample_val_valid(subsample_val):
    return (is_subsample_val_float(subsample_val) or is_subsample_val_int(subsample_val))

def create_non_exist_dir(dir_path):
    target_dir = Path(dir_path)
    if not target_dir.is_dir():
        target_dir.mkdir(mode=0o770, parents=True, exist_ok=True)

def extract_read_ids_for_target_chrom(align_file_path, target_chrom=None, min_map_qual=0, is_primary_align_only=False,
                                      is_include_unmapped_read=False, is_show_read_stats=True):
    align_file_ext = PurePosixPath(align_file_path).suffix
    if align_file_ext == '.bam':
        src_align = pysam.AlignmentFile(align_file_path, 'rb')
    elif align_file_ext == '.cram':
        src_align = pysam.AlignmentFile(align_file_path, 'rc')
    else:
        print(f'Alignment file type unknown for \'{align_file_path}\'.')
        return list()

    target_read_id_to_len_map = dict()
    read_count = 0

    for aligned_read in src_align:
        if aligned_read.is_unmapped:
            if not is_include_unmapped_read:
                continue
        else:
            if target_chrom is not None and aligned_read.reference_name != target_chrom:
                continue

            if is_primary_align_only and (aligned_read.is_secondary or aligned_read.is_supplementary):
                continue

            if aligned_read.mapping_quality < min_map_qual:
                continue

        read_len = aligned_read.infer_read_length()

        if aligned_read.query_name in target_read_id_to_len_map:
            if read_len > target_read_id_to_len_map[aligned_read.query_name]:
                target_read_id_to_len_map[aligned_read.query_name] = read_len
        else:
            target_read_id_to_len_map[aligned_read.query_name] = read_len

    if is_show_read_stats:
        if target_chrom is None:
            print(f'Number of target reads: {len(target_read_id_to_len_map.keys())}')
        else:
            print(f'Number of target reads for \'{target_chrom}\': {len(target_read_id_to_len_map.keys())}')

        #extract_read_lens = np.array(target_read_lens)
        extract_read_lens = np.array(list(target_read_id_to_len_map.values()))
        mean_read_len = round(np.mean(extract_read_lens), 2)
        sd_read_len = round(np.std(extract_read_lens), 2)
        max_read_len = np.max(extract_read_lens)
        min_read_len = np.min(extract_read_lens)

        print(f'Mean read length (S.D.) = {mean_read_len} ({sd_read_len})')
        print(f'Max read length = {max_read_len}')
        print(f'Min read length = {min_read_len}')

    src_align.close()

    target_read_ids = list(target_read_id_to_len_map.keys())
    target_read_ids.sort()

    return target_read_ids

def get_output_file_name_prefix(src_file_path, target_chrom=None, user_label=None):
    src_file_name_prefix = str(PurePosixPath(src_file_path).stem)

    if target_chrom is not None:
        src_file_name_prefix = f'{src_file_name_prefix}.{target_chrom}'

    if user_label is not None:
        src_file_name_prefix = f'{src_file_name_prefix}.{user_label}'

    return src_file_name_prefix

def generate_target_reads_file(target_read_ids, output_dir_path, output_reads_file_name_prefix):
    target_reads_file_path = str(PurePosixPath(output_dir_path).joinpath(f'{output_reads_file_name_prefix}.txt'))
    
    with open(target_reads_file_path, 'w') as f:
        for read_id in target_read_ids:
            f.write(f'{read_id}\n')

    return target_reads_file_path

def parse_target_reads_file(target_reads_file_path):
    with open(target_reads_file_path) as f:
        target_read_ids = list(f.readlines())

    target_read_ids = list(set(map(str.rstrip, target_read_ids)))
    target_read_ids.sort()

    return target_read_ids

def scan_dir_for_read_file_type(src_read_dir_paths):
    is_fast5_exist = False
    is_pod5_exist = False

    for dir_path in src_read_dir_paths:
        for child in Path(dir_path).iterdir():
            if child.is_dir():
                is_child_fast5_exist, is_child_pod5_exist = scan_dir_for_read_file_type([child])
                if is_child_fast5_exist and is_child_pod5_exist:
                    return is_child_fast5_exist, is_child_pod5_exist

                if is_child_fast5_exist:
                    is_fast5_exist = is_child_fast5_exist

                if is_child_pod5_exist:
                    is_pod5_exist = is_child_pod5_exist
            else:
                read_file_ext = PurePosixPath(child).suffix

                if read_file_ext == FAST5_FILE_EXT:
                    is_fast5_exist = True
                elif read_file_ext == POD5_FILE_EXT:
                    is_pod5_exist = True

            if is_fast5_exist and is_pod5_exist:
                return is_fast5_exist, is_pod5_exist

    return is_fast5_exist, is_pod5_exist

def subsample_reads(target_read_ids, subsample_val, seed=None):
    if seed is not None:
        random.seed(seed)

    target_read_count = len(target_read_ids)

    if subsample_val < 1:
        subsample_size = round(target_read_count * subsample_val)
    else:
        subsample_size = subsample_val

    if seed is None:
        print(f'Subsampling {subsample_size} reads out of {target_read_count} reads')
    else:
        print(f'Subsampling {subsample_size} reads out of {target_read_count} reads (seed: {seed})')

    return random.sample(target_read_ids, k=subsample_size)

def extract_target_reads_from_fast5_dir(fast5_dir_paths, output_dir_path, target_reads_file_path, output_reads_file_name_prefix, threads=1):
    for dir_path in fast5_dir_paths:
        if len(fast5_dir_paths) == 1:
            target_output_dir_path = output_dir_path
        else:
            target_output_dir_path = str(PurePosixPath(output_dir_path).joinpath(PurePosixPath(dir_path).stem))
            create_non_exist_dir(target_output_dir_path)

        cmd_str = f'python {ONT_FAST5_SUBSET_TOOL_PATH} -r -i {fast5_dir_paths} -s {target_output_dir_path} -l {target_reads_file_path} \
                    -f {output_reads_file_name_prefix}. -t {threads}'
        cmd_args = shlex.split(cmd_str)
        subprocess.run(cmd_args)

def extract_target_reads_from_pod5_dir(pod5_dir_paths, output_dir_path, target_reads_file_path, output_reads_file_name_prefix, threads=1):
    pod5_dir_path_str = ' '.join(pod5_dir_paths)
    output_pod5_file_path = str(PurePosixPath(output_dir_path).joinpath(f'{output_reads_file_name_prefix}.pod5'))
    cmd_str = f'pod5 filter -r {pod5_dir_path_str} -i {target_reads_file_path} -o {output_pod5_file_path} --missing-ok -t {threads}'
    cmd_args = shlex.split(cmd_str)
    subprocess.run(cmd_args)

def main():
    parser = ArgumentParser('Extraction tool for ONT long reads in FAST5/POD5 format')
    parser.add_argument('-i', '--input', required=True, action='store', dest='src_read_dir_paths', nargs='+',
                       help=f'Source long reads FAST5/POD5 directory paths')
    parser.add_argument('-o', '--output', required=True, action='store', dest='output_dir_path',
                        help='Output directory path for extracted long reads')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-b', '--bam', action='store', dest='bam_file_path', help='BAM file path')
    group.add_argument('-t', '--target', action='store', dest='read_id_file_path', help='Target read Ids file path')
    #parser.add_argument('-c', '--cram', required=False, help=f'CRAM file path (default: {DEFAULT_CRAM_FILE_PATH})',
    #                    default=DEFAULT_CRAM_FILE_PATH)
    parser.add_argument('-m', '--chrom', required=False, action='store', dest='chrom', help='Target chromosome (for BAM file only)')
    parser.add_argument('-p', '--primary', required=False, action='store_true', dest='is_primary_align_only',
                        help='Includes primary alignment only  (for BAM file only)', default=False)
    parser.add_argument('-u', '--unmapped', required=False, action='store_true', dest='is_include_unmapped_read',
                        help='Includes unmapped reads (overrides target chromosome, for BAM file only)', default=False)
    parser.add_argument('-l', '--label', required=False, action='store', dest='user_label',
                        help='User-defined label to be appended to the output reads file name')
    parser.add_argument('-s', '--subsample', required=False, action='store', dest='subsample_val',
                        help='Randomly subsample reads based on fraction (between 0 and 1) or absolute number (+ve integer > 1)')
    parser.add_argument('-S', '--seed', required=False, type=int, action='store', dest='seed',
                        help='Seed value for random subsampling')
    parser.add_argument('-d', '--dry-run', required=False, action='store_true', dest='dry_run',
                        help='Dry-run mode. Show reads statistics only')
    default_cpu_count = cpu_count()
    parser.add_argument('-T', '--threads', action='store', dest='threads', type=int, default=default_cpu_count,
                        help=f'Maximum no. of threads (default: {default_cpu_count})')
    args = parser.parse_args()

    for dir_path in args.src_read_dir_paths:
        if not Path(dir_path).is_dir():
            sys.exit(f'Source reads directory \'{dir_path}\' does not exist.')

    if args.chrom is not None and not is_input_chrom_valid(args.chrom):
        sys.exit(f'\'{args.chrom}\' is not a valid chromosome label.')

    if args.user_label is not None and not is_user_label_valid(args.user_label):
        sys.exit(f'\'{args.user_label}\' is not a valid user-defined file name label.')

    if args.subsample_val is None:
        subsample_val = 1
    else:
        if is_subsample_val_valid(args.subsample_val):
            if is_subsample_val_float(args.subsample_val):
                subsample_val = float(args.subsample_val)
            else:
                subsample_val = int(args.subsample_val)
        else:
            sys.exit('Subsampling value must be either a fraction between 0 and 1 (exclusive) or a positive integer.')

    is_fast5_exist, is_pod5_exist = scan_dir_for_read_file_type(args.src_read_dir_paths)
    if not is_fast5_exist and not is_pod5_exist:
        sys.exit(f'No reads files are found in \'{args.src_read_dir_paths}\'.')

    if args.threads < 1:
        proc_threads = default_cpu_count
    else:
        proc_threads = args.threads

    if not args.dry_run:
        '''
        if args.chrom is None:
            output_dir_path = args.output_dir_path
        else:
            output_dir_path = str(PurePosixPath(args.output_dir_path).joinpath(args.chrom))
        
        create_non_exist_dir(output_dir_path)
        '''

        create_non_exist_dir(args.output_dir_path)

    if args.bam_file_path is not None:
        target_read_ids = extract_read_ids_for_target_chrom(args.bam_file_path, target_chrom=args.chrom,
                                                            is_primary_align_only=args.is_primary_align_only,
                                                            is_include_unmapped_read=args.is_include_unmapped_read)
        if not args.dry_run:
            output_reads_file_name_prefix = get_output_file_name_prefix(args.bam_file_path, args.chrom, args.user_label)
            #target_reads_file_path = generate_target_reads_file(target_read_ids, output_dir_path, output_reads_file_name_prefix)
            target_reads_file_path = generate_target_reads_file(target_read_ids, args.output_dir_path, output_reads_file_name_prefix)
    elif args.read_id_file_path is not None:
        target_read_ids = parse_target_reads_file(args.read_id_file_path)
        target_read_count = len(target_read_ids)
        print(f'{target_read_count} unique reads are found in \'{args.read_id_file_path}\'')

        if target_read_count == 0 and not args.dry_run:
            sys.exit(f'Target read Id file \'{args.read_id_file_path}\' is empty.')

        target_reads_file_path = args.read_id_file_path

        if not args.dry_run:
            output_reads_file_name_prefix = get_output_file_name_prefix(args.read_id_file_path, args.chrom, args.user_label)
    else:
        sys.exit('No BAM file or target read Id file is provided.')

    if args.dry_run:
        print('No action performed in dry-run mode.')
        sys.exit(0)

    if subsample_val == 1:
        output_read_ids = target_read_ids
        final_target_read_file_path = target_reads_file_path
    else:
        output_read_ids = subsample_reads(target_read_ids, subsample_val, args.seed)
        subsample_size = len(output_read_ids)
        output_reads_file_name_prefix = f'{output_reads_file_name_prefix}.subsample_{subsample_size}'
        #final_target_read_file_path = generate_target_reads_file(output_read_ids, output_dir_path, output_reads_file_name_prefix)
        final_target_read_file_path = generate_target_reads_file(output_read_ids, args.output_dir_path, output_reads_file_name_prefix)

    print('Exporting reads...')

    if is_fast5_exist:
        #extract_target_reads_from_fast5_dir(args.src_read_dir_paths, output_dir_path, final_target_read_file_path,
        #                                    output_reads_file_name_prefix, proc_threads)
        extract_target_reads_from_fast5_dir(args.src_read_dir_paths, args.output_dir_path, final_target_read_file_path,
                                            output_reads_file_name_prefix, proc_threads)
    
    if is_pod5_exist:
        #extract_target_reads_from_pod5_dir(args.src_read_dir_paths, output_dir_path, final_target_read_file_path,
        #                                   output_reads_file_name_prefix, proc_threads)
        extract_target_reads_from_pod5_dir(args.src_read_dir_paths, args.output_dir_path, final_target_read_file_path,
                                           output_reads_file_name_prefix, proc_threads)

    #print(f'Extracted reads exported to directory \'{output_dir_path}\'')
    print(f'Extracted reads exported to directory \'{args.output_dir_path}\'')

    if args.bam_file_path is not None and not args.dry_run and final_target_read_file_path is not None:
        print(f'The read Id list is exported to {final_target_read_file_path}')

if __name__ == '__main__':
    main()
