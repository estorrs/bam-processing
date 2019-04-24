import argparse
import os
import subprocess

import bam_processing as bp

parser = argparse.ArgumentParser()

parser.add_argument('input_bam', type=str,
        help='input bam')

parser.add_argument('--reference-fasta', type=str,
        help='reference fasta')
parser.add_argument('--known-sites', type=str,
        help='known sites for base recal')
parser.add_argument('--output', type=str,
        default='output.bam', help='output bam')
parser.add_argument('--workflow-type', type=str,
        default='standard', help='Which workflow to run')
parser.add_argument('--fixmate', action='store_true',
        help='run samtools fixmate as part of workflow.')
parser.add_argument('--properly-paired-only', action='store_true',
        help='filter out all alignments that are not properly paired.')
parser.add_argument('--fix-255-mapping-quality', action='store_true',
        help='Changes all 255 mapping qualities to 60.')
parser.add_argument('--temp-files-dir', type=str,
        help='directory to put samtools temporary files in.')

args = parser.parse_args()

def run_standard_workflow(input_bam, output_bam, reference_fasta, known_sites_vcf_gz,
        fixmate, properly_paired_only, fix_255_mapping_quality, temp_files_dir):
    if temp_files_dir is None:
        bp.run_basic_preprocessing(input_bam, output_bam, reference_fasta, known_sites_vcf_gz,
                fixmates=fixmate, properly_paired_only=properly_paired_only)
    else:
        bp.run_basic_preprocessing(input_bam, output_bam, reference_fasta, known_sites_vcf_gz,
                fixmates=fixmate, properly_paired_only=properly_paired_only,
                fix_255_mapping_quality=fix_255_mapping_quality, temp_files_dir=temp_files_dir)

def main():
    run_standard_workflow(args.input_bam, args.output, args.reference_fasta, args.known_sites,
            args.fixmate, args.properly_paired_only, args.fix_255_mapping_quality,
            args.temp_files_dir)

if __name__ == '__main__':
    main()
