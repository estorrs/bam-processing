import os
import re
import shutil
import subprocess
import uuid

def create_sorted_bam(bam_fp, output_fp=None, name_sorted=False, max_memory='10G',
        temp_files_dir=os.getcwd(), use_temp_for_output=True):
    """Creates sorted bam and returns the filepath.

    If no filepath given then default is used.

    If bam is already sorted then original filepath is returned

    Returns
        fp - filepath for sorted bam
    """
    tool_args = ('samtools', 'stats', bam_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')
 
    # if already sorted then return input filepath
    # ignore this if trying to name sort
    if 'is sorted:\t1' in output and not name_sorted:
        return bam_fp

    if output_fp is None:
        if use_temp_for_output:
            output_fp = os.path.join(temp_files_dir, f'sorted.{str(uuid.uuid4())}.bam')
        else:
            output_fp = f'sorted.{str(uuid.uuid4())}.bam'

    if name_sorted:
        tool_args = ('samtools', 'sort', '-n', '-m', max_memory, '-T',
                os.path.join(temp_files_dir, str(uuid.uuid4())), '-o', output_fp, bam_fp)
    else:
        tool_args = ('samtools', 'sort', '-m', max_memory, '-T',
                os.path.join(temp_files_dir, str(uuid.uuid4())), '-o', output_fp, bam_fp)

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

def run_add_or_replace_read_groups(input_fp, output_fp, temp_files_dir=os.getcwd()):
    # sort and index if needed
    sorted_input_fp = create_sorted_bam(input_fp, temp_files_dir=temp_files_dir)
    index_bam(sorted_input_fp)

    tool_args = add_or_replace_read_groups(input_fp=sorted_input_fp, output_fp=output_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')

    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)
        os.remove(sorted_input_fp + '.bai')

    return output

def split_n_cigar_reads(reference_fp, input_fp='/dev/stdin', output_fp='/dev/stdout'):
    tool_args = ('gatk', 'SplitNCigarReads',
            '-R', reference_fp,
            '-I', input_fp,
            '-O', output_fp)

    return tool_args

def run_split_n_cigar_reads(input_fp, output_fp, reference_fp,
        temp_files_dir=os.getcwd()):
    # sort and index if needed
    sorted_input_fp = create_sorted_bam(input_fp, temp_files_dir=temp_files_dir)
    index_bam(sorted_input_fp)

    # make sure reference is prepared
    index_reference(reference_fp)

    tool_args = split_n_cigar_reads(reference_fp, input_fp=sorted_input_fp, output_fp=output_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')

    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)
        os.remove(sorted_input_fp + '.bai')

    return output

def mark_duplicates(input_fp='/dev/stdin', output_fp='/dev/stdout', max_records_in_ram=100000,
        temp_dir='temp_dir', metrics_fp='output.metrics'):
    tool_args = ('picard', 'MarkDuplicates',
            f'I={input_fp}',
            f'O={output_fp}',
            f'MAX_RECORDS_IN_RAM={max_records_in_ram}',
            'VALIDATION_STRINGENCY=SILENT',
            f'TMP_DIR={temp_dir}',
            f'M={metrics_fp}')

    return tool_args

def run_mark_duplicates(input_fp, output_fp, temp_files_dir=os.getcwd()):
    # sort and index if needed
    sorted_input_fp = create_sorted_bam(input_fp, temp_files_dir=temp_files_dir)
    index_bam(sorted_input_fp)
    metrics_fp = os.path.join(temp_files_dir, f'output.{uuid.uuid4()}.metrics')
    temp_dir = os.path.join(temp_files_dir, f'temp_{uuid.uuid4()}_dir')
    os.mkdir(temp_dir)

    tool_args = mark_duplicates(input_fp=sorted_input_fp, output_fp=output_fp,
            temp_dir=temp_dir, metrics_fp=metrics_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')

    os.remove(metrics_fp)
    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)
        os.remove(sorted_input_fp + '.bai')

    # remove temporary directory
    shutil.rmtree(temp_dir)

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

def run_base_recalibration(input_fp, output_fp, reference_fp, known_sites_fp,
        temp_files_dir=os.getcwd()):
    # make sure reference is prepared
    index_reference(reference_fp)
    create_reference_sequence_dict(reference_fp)

    # sort and index bam if needed
    sorted_input_fp = create_sorted_bam(input_fp, temp_files_dir=temp_files_dir)
    index_bam(sorted_input_fp)

    table_fp = os.path.join(temp_files_dir, f'output.{str(uuid.uuid4())}.table')
    tool_args = base_recalibrator_table(sorted_input_fp, table_fp, reference_fp, known_sites_fp)
    output = subprocess.check_output(tool_args).decode('utf-8')

    tool_args = base_recalibration(sorted_input_fp, output_fp, reference_fp, table_fp)
    output += '\n\n' + subprocess.check_output(tool_args).decode('utf-8')

    os.remove(table_fp)
    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)
        os.remove(sorted_input_fp + '.bai')

    return output

def fixmates(input_fp='/dev/stdin', output_fp='/dev/stdout'):
    tool_args = ('samtools', 'fixmate', input_fp, output_fp)

    return tool_args

def run_fixmates(input_fp, output_fp, temp_files_dir=os.getcwd()):
    # sort and index bam if needed
    sorted_input_fp = create_sorted_bam(input_fp, name_sorted=True, temp_files_dir=temp_files_dir)

    tool_args = fixmates(input_fp=sorted_input_fp, output_fp=output_fp)
    subprocess.check_output(tool_args)

    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)

def properly_paired(input_fp='/dev/stdin', output_fp='/dev/stdout'):
    tool_args = ('samtools', 'view', '-h', '-f', '0x2', '-o', output_fp, input_fp)

    return tool_args

def run_properly_paired(input_fp, output_fp, temp_files_dir=os.getcwd()):
    # sort and index bam if needed
    sorted_input_fp = create_sorted_bam(input_fp, temp_files_dir=temp_files_dir)
    index_bam(sorted_input_fp)

    tool_args = properly_paired(input_fp=sorted_input_fp, output_fp=output_fp)
    subprocess.check_output(tool_args)

    if sorted_input_fp != input_fp:
        os.remove(sorted_input_fp)
        os.remove(sorted_input_fp + '.bai')

def run_fix_255_mapping_quality(input_fp, output_fp):
    ps_1 = subprocess.Popen(('samtools', 'view', '-h', input_fp), stdout=subprocess.PIPE)
    ps_2 = subprocess.Popen(('sed', r's/^\([^	]*	[^	]*	[^	]*	[^	]*	\)255/\160/'), stdin=ps_1.stdout, stdout=subprocess.PIPE)
    ps_3 = subprocess.Popen(('sed', r's/MQ:i:255/MQ:i:60/'), stdin=ps_2.stdout, stdout=subprocess.PIPE)

    subprocess.check_output(('samtools', 'view', '-h', '-o', output_fp), stdin=ps_3.stdout)
    ps_1.wait()
    ps_2.wait()
    ps_3.wait()

def run_cptac3_preprocessing(input_fp, output_fp, reference_fp, temp_files_dir=os.getcwd()):
    # sort and add readgroups
    read_group_output = os.path.join(temp_files_dir, f'read_groups.{str(uuid.uuid4())}.bam')
    run_add_or_replace_read_groups(input_fp, read_group_output,
            temp_files_dir=temp_files_dir)

    # mark duplicates
    mark_duplicates_output = os.path.join(temp_files_dir, f'mark_duplicates.{str(uuid.uuid4())}.bam')
    run_mark_duplicates(read_group_output, mark_duplicates_output,
            temp_files_dir=temp_files_dir)
    # remove temp output
    os.remove(read_group_output)
    os.remove(read_group_output + '.bai')

    # split n cigar reads
    split_output = os.path.join(temp_files_dir, f'split.{str(uuid.uuid4())}.bam')
    run_split_n_cigar_reads(mark_duplicates_output, output_fp, reference_fp,
            temp_files_dir=temp_files_dir)
    # remove temp output
    os.remove(mark_duplicates_output)
    os.remove(mark_duplicates_output + '.bai')
    os.remove(output_fp.replace('.bam', '.bai'))

    sorted_output = create_sorted_bam(output_fp, temp_files_dir=temp_files_dir)
    index_bam(sorted_output)
    shutil.move(sorted_output, output_fp)

def run_basic_preprocessing(input_fp, output_fp, reference_fp, known_sites_fp,
        properly_paired_only=False, fixmates=False, temp_files_dir=os.getcwd(),
        fix_255_mapping_quality=False):
    original_input_fp = input_fp

    if fixmates:
        fixmates_output = os.path.join(temp_files_dir, f'temp.{str(uuid.uuid4())}.bam')
        run_fixmates(input_fp, fixmates_output, temp_files_dir=temp_files_dir)
        input_fp = fixmates_output

    if properly_paired_only:
        properly_paired_output = os.path.join(temp_files_dir, f'temp.{str(uuid.uuid4())}.bam')
        run_properly_paired(input_fp, properly_paired_output, temp_files_dir=temp_files_dir)
        if input_fp != original_input_fp:
            os.remove(input_fp)
            if os.path.isfile(input_fp + '.bai'):
                os.remove(input_fp + '.bai')
        input_fp = properly_paired_output

    if fix_255_mapping_quality:
        fix_255_output = os.path.join(temp_files_dir, f'temp.{str(uuid.uuid4())}.bam')
        run_fix_255_mapping_quality(input_fp, fix_255_output)
        if input_fp != original_input_fp:
            os.remove(input_fp)
            if os.path.isfile(input_fp + '.bai'):
                os.remove(input_fp + '.bai')
        input_fp = fix_255_output

    # add read groups
    read_group_output = os.path.join(temp_files_dir, f'temp.{str(uuid.uuid4())}.bam')
    run_add_or_replace_read_groups(input_fp, read_group_output, temp_files_dir=temp_files_dir)
    if input_fp != original_input_fp:
        os.remove(input_fp)
        if os.path.isfile(input_fp + '.bai'):
            os.remove(input_fp + '.bai')

    # mark duplicates
    mark_duplicates_output = os.path.join(temp_files_dir, f'temp.{str(uuid.uuid4())}.bam')
    run_mark_duplicates(read_group_output, mark_duplicates_output,
            temp_files_dir=temp_files_dir)
    # remove temp output
    os.remove(read_group_output)
    os.remove(read_group_output + '.bai')

    # base recalibration
    run_base_recalibration(mark_duplicates_output, output_fp, reference_fp, known_sites_fp,
            temp_files_dir=temp_files_dir)
    # remove temp output
    os.remove(mark_duplicates_output)
    os.remove(mark_duplicates_output + '.bai')
    os.remove(output_fp.replace('.bam', '.bai'))

    sorted_output = create_sorted_bam(output_fp, temp_files_dir=temp_files_dir)
    index_bam(sorted_output)

    shutil.move(sorted_output, output_fp)
