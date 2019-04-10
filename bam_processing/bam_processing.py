import os
import re
import shutil
import subprocess
import uuid

def create_sorted_bam(bam_fp, output_fp=None):
    """Creates sorted bam and returns the filepath.

    If no filepath given then default is used.

    If bam is already sorted then original filepath is returned

    Returns
        fp - filepath for sorted bam
    """
    tool_args = ('samtools', 'stats', bam_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')
    
    # if already sorted then return input filepath
    if 'is sorted:\t1' in output:
        return bam_fp

    if output_fp is None:
        output_fp = f'sorted.{str(uuid.uuid4())}.bam'
    tool_args = ('samtools', 'sort', '-o', output_fp, bam_fp)
    subprocess.check_output(tool_args)

    return output_fp

def index_bam(bam_fp):
    """index the given bam if it is not already"""
    if not os.path.isfile(f'{bam_fp}.bai'):
        tool_args = ('samtools', 'index', bam_fp)
        subprocess.check_output(tool_args)

def index_reference(reference_fp):
    if not os.path.isfile(f'{reference_fp}.fai'):
        tool_args = ('samtools', 'faidx', reference_fp)
        subprocess.check_output(tool_args)

def create_reference_sequence_dict(reference_fp):
    index_reference(reference_fp)

    output_fp = re.sub(r'.[^.]*$', '.dict', reference_fp)
    if not os.path.isfile(output_fp):
        tool_args = ('picard', 'CreateSequenceDictionary',
                f'R={reference_fp}',
                'O=' + output_fp,
                )
        subprocess.check_output(tool_args)

def index_vcf(bgzip_vcf_fp):
    if not os.path.isfile(f'{bgzip_vcf_fp}.tbi'):
        tool_args = ('tabix', '-p', 'vcf', bgzip_vcf_fp)
        subprocess.check_output(tool_args)


def add_or_replace_read_groups(input_fp='/dev/stdin', output_fp='/dev/stdout'):
    tool_args = ('picard', 'AddOrReplaceReadGroups',
            f'I={input_fp}',
            f'O={output_fp}',
            'SO=coordinate', 'RGID=id', 'RGLB=library',
            'RGPL=platform', 'RGPU=machine', 'RGSM=sample')

    return tool_args

def run_add_or_replace_read_groups(input_fp, output_fp):
    # sort and index if needed
    sorted_input_fp = create_sorted_bam(input_fp)
    index_bam(sorted_input_fp)

    tool_args = add_or_replace_read_groups(input_fp=sorted_input_fp, output_fp=output_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')

    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)

    return output

def mark_duplicates(input_fp='/dev/stdin', output_fp='/dev/stdout', max_records_in_ram=100000,
        metrics_fp='output.metrics'):
    tool_args = ('picard', 'MarkDuplicates',
            f'I={input_fp}',
            f'O={output_fp}',
            f'MAX_RECORDS_IN_RAM={max_records_in_ram}',
            'VALIDATION_STRINGENCY=SILENT',
            f'M={metrics_fp}')

    return tool_args

def run_mark_duplicates(input_fp, output_fp):
    # sort and index if needed
    sorted_input_fp = create_sorted_bam(input_fp)
    index_bam(sorted_input_fp)
    metrics_fp = f'output.{uuid.uuid4()}.metrics'

    tool_args = mark_duplicates(input_fp=sorted_input_fp, output_fp=output_fp,
            metrics_fp=metrics_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')

    os.remove(metrics_fp)
    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)

    return output

def base_recalibrator_table(input_fp, output_fp, reference_fp, known_sites_fp):
    tool_args = ('gatk', 'BaseRecalibrator',
            '-I', input_fp,
            '-R', reference_fp,
            '--known-sites', known_sites_fp,
            '-O', output_fp)

    return tool_args

def base_recalibration(input_fp, output_fp, reference_fp, table_fp):
    tool_args = ['gatk', 'ApplyBQSR',
            '-I', input_fp,
            '-R', reference_fp,
            '--bqsr-recal-file', table_fp,
            '-O', output_fp]
    return tool_args

def run_base_recalibration(input_fp, output_fp, reference_fp, known_sites_fp):
    # make sure reference is prepared
    index_reference(reference_fp)
    create_reference_sequence_dict(reference_fp)

#     # make sure known sites is indexed
#     index_vcf(known_sites_fp)

    # sort and index bam if needed
    sorted_input_fp = create_sorted_bam(input_fp)
    index_bam(sorted_input_fp)

    table_fp = f'output.{str(uuid.uuid4())}.table'
    tool_args = base_recalibrator_table(sorted_input_fp, table_fp, reference_fp, known_sites_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')

    tool_args = base_recalibration(sorted_input_fp, output_fp, reference_fp, table_fp)
    output += '\n\n' + subprocess.check_output(tool_args).decode('utf-8')

    os.remove(table_fp)
    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)

    return output

def run_basic_preprocessing(input_fp, output_fp, reference_fp, known_sites_fp):
    # add read groups
    read_group_output = f'temp.{str(uuid.uuid4())}.bam'
    run_add_or_replace_read_groups(input_fp, read_group_output)

    # mark duplicates
    mark_duplicates_output = f'temp.{str(uuid.uuid4())}.bam'
    run_mark_duplicates(read_group_output, mark_duplicates_output)
    # remove temp output
    os.remove(read_group_output)
    os.remove(read_group_output + '.bai')

    # base recalibration
    run_base_recalibration(mark_duplicates_output, output_fp, reference_fp, known_sites_fp)
    # remove temp output
    os.remove(mark_duplicates_output)
    os.remove(mark_duplicates_output + '.bai')
    os.remove(output_fp.replace('.bam', '.bai'))

    sorted_output = create_sorted_bam(output_fp)
    index_bam(sorted_output)

    shutil.move(sorted_output, output_fp)

# def base_recalibration(input_fp, reference_fp, known_sites_fp, output_fp='/dev/stdout'):
#     # prepare reference and known sites

