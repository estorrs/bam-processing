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

args = parser.parse_args()

def run_standard_workflow(input_bam, output_bam, reference_fasta, known_sites_vcf_gz):
    bp.run_basic_preprocessing(input_bam, output_bam, reference_fasta, known_sites_vcf_gz)

def main():
    run_standard_workflow(args.input_bam, args.output, args.reference_fasta, args.known_sites)

if __name__ == '__main__':
    main()
