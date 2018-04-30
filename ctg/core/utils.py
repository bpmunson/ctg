"""
CTG utility functions.
"""

import numpy as np
import gzip
import pysam
from collections import defaultdict

def readfq(fp): # this is a generator function
    """ Fastq file reader from Heng Li: https://github.com/lh3/readfq/blob/master/readfq.py 

        Args: 
            fp (file handle): a file handle to the fastq file
        Yields:
            read (tuple):
                (read_id (str), sequence (str), quality (str) ) 
                if the read is a fasta record then no quality is included in the read tuple

    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def consensus(fastq, threshold=0.75, num_reads=1000):
    """ Calculate a quick consensus sequence from a fastq file.

        Args:
            fastq (str): path to fastq file to take reads from
            threshold (float): the frequency threshold to call a base
            num_reads (int): use the first n reads of the fastq file
        Returns: None
    """
    c = 0
    x = []

    if fastq.endswith("gz"):
        with gzip.open(fastq, "rb") as handle:
            fq_reader = readfq(handle)
            while c<num_reads:
                x.append(next(fq_reader)[1])
                c+=1
    else: 
        with open(fastq, "r") as handle:
            fq_reader = readfq(handle)
            while c<num_reads:
                x.append(next(fq_reader)[1])
                c+=1

    d = defaultdict(lambda: defaultdict(int))
    longest = max([len(s) for s in x])
    for i in range(longest):
        for s in x:
            d[i][s[i]]+=1
    consensus = []
    hs = []
    # calculate frequencies
    for i in d:
        called = False
        tot = sum(d[i].values())
        h = 0
        for j in d[i]:
            f = d[i][j]/tot
            #h += -f*np.log2(f) 
            if f>threshold:
                consensus.append(j)
                called = True
                break
        if not called:
            consensus.append("N")
    consensus = "".join(consensus)
    print(consensus)

def mate_pair_bam_reader(bam_file_path, paired=True):
    """
    Read mate pairs from a bam file.

    Args:
        bam_file_path (str): path to bam file to read
        paired (bool): True if the bam file contains paried reads, false
            otherwise
    Yeilds:
        read pair (pysam AlignedSegement, pysamAlignedSegment): a read pair,
            first is read 1, second is read 2 

    """
    if bam_file_path.endswith('bam'):
        method = "rb"
    else:
        method = "r"

    bam = pysam.AlignmentFile(bam_file_path)

    while True:
        # get a read
        try:
            read = next(bam)
        except StopIteration:
            break

        # get pointer
        pointer = bam.tell()

        if paired:
            # if we expect paired reads, try and get the mate, ignoring unpaired reads for now
            try:
                mate = next(bam)
            except StopIteration:
                break #
            if read.query_name != mate.query_name:
                # these are not a pair
                # so return to previous position and try again
                bam.seek(pointer)
                continue
        else:
            mate = None
        yield read, mate

def get_fragment_count(bam):
    """ 
    Get the number of fragments an aligned bam, checking for paired status
    first.

    Args: 
        bam (pysam.AlignmentFile): sam/bam file to get fragment count of
    Returns:
        readcount (int): the number of fragments in the bam file
        True/False bool:  True if the fragments are paired, False otherwise
    """
    # check paired
    read_count = int(pysam.view("-c", "-f65", bam).rstrip())
    if read_count == 0:
        # get singles
        read_count = int(pysam.view("-c",bam).rstrip())
        return read_count, False
    else:
        return read_count, True
