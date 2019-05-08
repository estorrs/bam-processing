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
parser.add_argument('--max-memory', type=str,
        default='1g', help='max heap size for java to allocate')

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

def run_cptac3_workflow(input_bam, output_bam, reference_fasta, temp_files_dir, max_memory):
    if temp_files_dir is None:
        bp.run_cptac3_preprocessing(input_bam, output_bam, reference_fasta, max_mem=max_memory)
    else:
        bp.run_cptac3_preprocessing(input_bam, output_bam, reference_fasta, max_mem=max_memory,
                temp_files_dir=temp_files_dir)

def run_cptac2_workflow(input_bam, output_bam, reference_fasta, temp_files_dir, max_memory):
    if temp_files_dir is None:
        bp.run_cptac2_preprocessing(input_bam, output_bam, reference_fasta, max_mem=max_memory)
    else:
        bp.run_cptac2_preprocessing(input_bam, output_bam, reference_fasta, max_mem=max_memory,
                temp_files_dir=temp_files_dir)

def main():
    if args.workflow_type == 'standard':
        run_standard_workflow(args.input_bam, args.output, args.reference_fasta, args.known_sites,
                args.fixmate, args.properly_paired_only, args.fix_255_mapping_quality,
                args.temp_files_dir)
    elif args.workflow_type == 'cptac3':
        run_cptac3_workflow(args.input_bam, args.output, args.reference_fasta, args.temp_files_dir,
                args.max_memory)
    elif args.workflow_type == 'cptac2':
        run_cptac2_workflow(args.input_bam, args.output, args.reference_fasta, args.temp_files_dir,
                args.max_memory)
    else:
        raise ValueError('must specify correct workflow')

if __name__ == '__main__':
    main()
