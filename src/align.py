"""
Wrapper for bowtie2 functions
"""
import os
import subprocess
import logging

# Command line functions
def bowtie2_build(reference, bt2_index_base, *args, **kwargs):
    """ Function to call bowtie2-build from python. 

    Providing a valid reference fasta file, use subprocess to build the index 
    files necessary for bowtie2 alignment.

    No validation of arguments is performed.
    Args:
        reference (str): Path to a reference fasta file to build 
            bowtie2 index from
        bt2_index_base (str): Base path to write the index files to
        *args: additional bowtie2 arguments, examples: 
        **kwargs: additional bowtie2-build arguments, 
            examples: o=5, ftabchars=10
    Returns:
        None
    Raises:
        OSError is bowtie2 is not installed.
    """

    try:
        subprocess.check_output('bowtie2-build --version', shell=True)
    except OSError:
        raise OSError("{} {}".format(
            'Commnad bowtie2-build not found.',
            'Please install bowtie2.'))

    call = ["bowtie2-build"]
    # add optional arguments
    for arg in args:
        call.append(arg)
    for k, v in kwargs.items():
        if len(k)==1: # single 
            call.append("-{} {}".format(k, v))
        else:
            call.append("--{} {}".format(k, v))
    if bt2_index_base:
        call.append(bt2_index_base)
    call.append(reference_in)
    #print("call: ", call)
    logging.info("Building Bowtie2 Index.")
    subprocess.check_call(call)


    return 0

def bowtie2_align(fastq, bt2_index_base, *args, 
    guide_start=None, guide_length=None,
    expected_barcode=None, barcode_start=None,
    binary_flags_to_add = None, bam=None, **kwargs):
    """ 
    Align single reads containing a guide sequence with bowtie2 the guide
    based reference.  Additionally extracting the guide sequence,
    calculating the edit distance in this guide and add both as tags to bam 
    file containing the alignments.

    Optioanlly, extract a barcode sequence and calculate the edit distance, 
    again adding both as bam tags.

    Optionally, add bitwise flags to the alingments.

    Args:
        reference (str): Path to a reference fasta file to build bowtie2 
            index from
        bt2_index_base (str): Base path containing the bowtie2 index
        guide_start (int): Position from the 5' end of the read the guide
            sequence starts in the reads 
        guide_length (int): The expected length of the guide sequence
        expected_barcode: The expected barcode sequence in IUPAC format. 
            ex: WSWS for [A|T][G|C][A|T][G|C][A|T][G|C]
                so AGTC would be extracted
        barcode_start (int): Position from the 5' end of the read the barcode
            sequence starts in the reads 
        binary_flags_to_add (int): Bitwise flags to add to the reads. 
            ex: 65 for a paired read 1
        bam: the output bam file to write the single end reads to
        *args: additional bowtie2 arguments, examples: 
        **kwargs: additional bowtie2-build arguments,
             examples: o=5, ftabchars=10
    Returns:
        None
    Raises:
        OSError is bowtie2 is not installed.
        OSError is samtools is not installed.
    """
    try:
        subprocess.check_output('bowtie2 --version', shell=True)
    except OSError:
        raise OSError('Commnad bowtie2 not found. Please install bowtie2.')

    try:
        subprocess.check_output('samtools --version', shell=True)
    except OSError:
        raise OSError('Commnad samtools not found. Please install santools.')
    
    call = ['bowtie2']
    # add optional arguments
    for arg in args:
        call.append("{}".format(arg))
    for k, v in kwargs.items():
        call.append("--{} {}".format(k, v))

    call.append("-x {}".format(bt2_index_base))
    call.append("-U {}".format(fastq))

    # add 
    if guide_start or expected_barcode:
        call.append("| add_tags.py")
        if guide_start:
            if guide_length is None:
                raise RuntimeError()
            call.append("--guide_start {} --guide_length {}".format(
                guide_start,
                guide_length))
        if expected_barcode:
            if barcode_start is None:
                raise RuntimeError()
            call.append("--expected_barcode {} --barcode_start {}".format(
                expected_barcode,
                barcode_start))
        if binary_flags_to_add:
            call.append("--flag {}".format(binary_flags_to_add))


    # add conversion to bam format
    if bam is None:
        bam = "{}.bam".format(os.path.splitext(fastq)[0])

    call.append("| samtools view -Sb - > {}".format(bam))

    # merge call and execute    
    call = " ".join(call)
    logging.info("Running Bowtie2 Alignment: {}".format(call))
    subprocess.check_output(call, shell=True)

    return 0

def merge_fixmate(r1_bam, r2_bam, out_bam, tmp="./sort"):
    """ 
    Combine two single end alignments, sorting, and fixing mate information
    in the process.

    Uses subprocess to call samtools command line functions.
    
    Args:
        a1_bam (str): Input path to the read1 aligned bam 
        r2_bam (str): Input path to the read2 aligned bam 
        out_bam (str): Output path to the combined bam to
        tmp_dir (str): Temporary working directory for samtools to use
    Returns:
        None
    Raises:
        OSError is samtools is not installed.
    """
    try:
        subprocess.check_output('samtools --version', shell=True)
    except OSError:
        raise OSError('Commnad bowtie2 not found. Please install samtools.')   
    call = [    "samtools cat", r1_bam, r2_bam,
                "| samtools sort", "-T", tmp, "-n -",
                " | samtools fixmate -", out_bam]

    # merge call and execute    
    call = " ".join(call)
    logging.info("Merging and Fixing Mate Information: {}".format(call))
    subprocess.check_output(call, shell=True)
    return 0
