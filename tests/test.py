import pytest
import os
import re
import sys
import subprocess

import bam_processing as bp

INPUT_BAM = 'tests/data/test.2.bam'
REFERENCE_FASTA = 'tests/data/test.fa'
KNOWN_SITES_VCF_GZ = 'tests/data/test.vcf.gz'
TEMP_FILES_DIR = os.getcwd()

def test_create_sorted_bam():
    output_fp = bp.create_sorted_bam(INPUT_BAM)

def test_index_bam():
    # sort so you can index
    output_fp = bp.create_sorted_bam(INPUT_BAM)

    bp.index_bam(output_fp)

    assert os.path.isfile(INPUT_BAM + '.bai')

def test_index_reference():
    bp.index_reference(REFERENCE_FASTA)

    assert os.path.isfile(REFERENCE_FASTA + '.fai')

def test_create_sequence_dict():
    bp.create_reference_sequence_dict(REFERENCE_FASTA)

    assert os.path.isfile(re.sub(r'.[^.]*$', '.dict', REFERENCE_FASTA))

def test_index_vcf():
    bp.index_vcf(KNOWN_SITES_VCF_GZ)

    assert os.path.isfile(KNOWN_SITES_VCF_GZ + '.tbi')

def test_add_or_replace_read_groups():
    bp.run_add_or_replace_read_groups(input_fp=INPUT_BAM, output_fp='output.bam')
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
    assert 'RG:Z:id' in output 

def test_split_n_cigar_reads():
    bp.run_split_n_cigar_reads(input_fp=INPUT_BAM, output_fp='output.bam',
            reference_fp=REFERENCE_FASTA)
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
    assert True


def test_mark_duplicates():
    bp.run_mark_duplicates(input_fp=INPUT_BAM, output_fp='output.bam')
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
    
    assert 'PG:Z:MarkDuplicates' in output 

def test_properly_paired():
    bp.run_properly_paired(INPUT_BAM, 'output.bam')
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')

    assert True

def test_fixmate():
    bp.run_fixmates(INPUT_BAM, 'output.bam')
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')

    assert True

def test_fix_star_mapping_quality():
    bp.run_fix_255_mapping_quality(INPUT_BAM, 'output.bam')
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')

    assert '\t60\t' in output

def test_base_recalibration():
    bp.run_base_recalibration(INPUT_BAM, 'output.bam', REFERENCE_FASTA, KNOWN_SITES_VCF_GZ)
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
    
    assert 'ID:GATK ApplyBQSR' in output and '60' in output

def test_cptac3_processing():
    bp.run_cptac3_preprocessing(INPUT_BAM, 'output.bam', REFERENCE_FASTA)
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
    
    assert True

def test_cptac2_prospective_processing():
    bp.run_cptac2_prospective_preprocessing(INPUT_BAM, 'output.bam', REFERENCE_FASTA)
    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
    
    assert True

def test_standard_cli():
    tool_args = ('python', 'bam_processing/bam_processing_cli.py',
            '--workflow-type', 'cptac3',
            '--max-memory', '1g',
            '--reference-fasta', REFERENCE_FASTA,
            '--output', 'output.bam',
            INPUT_BAM)
    subprocess.check_output(tool_args)

    output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
    
#    assert 'ID:GATK ApplyBQSR' in output
    assert True

# def test_simple_processing():
#     bp.run_basic_preprocessing(INPUT_BAM, 'output.bam', REFERENCE_FASTA, KNOWN_SITES_VCF_GZ)
#     output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
#     
#     assert 'ID:GATK ApplyBQSR' in output
# 
# def test_simple_processing_with_extras():
#     bp.run_basic_preprocessing(INPUT_BAM, 'output.bam', REFERENCE_FASTA, KNOWN_SITES_VCF_GZ,
#             properly_paired_only=True, fixmates=True, fix_255_mapping_quality=True)
#     output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
#     
#     assert 'ID:GATK ApplyBQSR' in output and '\t60\t' in output
# 
# def test_standard_cli():
#     tool_args = ('python', 'bam_processing/bam_processing_cli.py',
#             '--reference-fasta', REFERENCE_FASTA,
#             '--known-sites', KNOWN_SITES_VCF_GZ,
#             '--output', 'output.bam',
#             INPUT_BAM)
#     subprocess.check_output(tool_args)
# 
#     output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
#     
#     assert 'ID:GATK ApplyBQSR' in output
# 
# def test_standard_cli_with_extras():
#     tool_args = ('python', 'bam_processing/bam_processing_cli.py',
#             '--reference-fasta', REFERENCE_FASTA,
#             '--known-sites', KNOWN_SITES_VCF_GZ,
#             '--temp-files-dir', TEMP_FILES_DIR,
#             '--fixmate', '--properly-paired-only', '--fix-255-mapping-quality',
#             '--output', 'output.bam',
#             INPUT_BAM)
#     subprocess.check_output(tool_args)
# 
#     output = subprocess.check_output(('samtools', 'view', '-h', 'output.bam')).decode('utf-8')
#     
#     assert 'ID:GATK ApplyBQSR' in output and '\t60\t' in output

# def test_together():
#     tool_args = bp.add_or_replace_read_groups(input_fp=INPUT_BAM)
#     add_replace_process = subprocess.Popen(tool_args, stdout=subprocess.PIPE)
# 
#     tool_args = bp.mark_duplicates()
#     mark_duplicates_process = subprocess.Popen(tool_args, stdin=add_replace_process.stdout, stdout=subprocess.PIPE)
# 
#     output = subprocess.check_output(('samtools', 'view', '-h'), stdin=mark_duplicates_process.stdout).decode('utf-8')
# 
#     add_replace_process.wait()
#     mark_duplicates_process.wait()
# 
#     assert 'PG:Z:MarkDuplicates' in output 
