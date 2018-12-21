"""
Functions to count constructs from fastq files.
"""

import os
import pysam
import logging
#import pandas as pd
import numpy as np
from collections import defaultdict
from ctg.core.utils import mate_pair_bam_reader
from ctg.core.utils import get_fragment_count



_construct_divider = "__"
_probe_divider = "-"

# get complement of a sequence
comp = lambda x: x.upper().translate(str.maketrans('ACGT','TGCA'))

# get reverse complement of seqeunce
rev_comp = lambda x: comp(x)[::-1]

def get_construct_divider():
    return _construct_divider

def get_probe_divider():
    return _probe_divider

def read_1_callback(read):
    """ Function to determine if read is the first in a pair
    
        Args:
            pysam AlignedSegment
        Return:
            True if bitwise flag 1 is set for read, False otherwise
    """
    if read.is_read1:
        return True
    else:
        return False

def read_library_file(library, sep="\t", header=0, comment="#"):
    """ Import a library file 
    Expected format is:
    construct_id    target_a_id probe_a_id  probe_a_sequence  ..
    target_b_id    probe_b_id  probe_b_sequence

    constructs_id is 
    Args:
        libary (str): path to library defintion
        sep (str): field delimiter
        header (str): header line after initial comments
        comment (str): character at begining of line specifying a comment
    Returns: 
        ld (dict): dictionary wtih construct_ids as keys and dictionary values
            dictionary values are of the form: column_id: row_value
    Raises:
        RuntimeError if duplicate constructs were seen in library 
    """
    ld = dict()
    options = dict()
    with open(library, 'r') as handle:
        c=0
        for line in handle:
            if line.startswith(comment):
                # expect tags
                try:
                    tag, value = line.lstrip("#").split(":")
                    tag = tag.lstrip()
                except ValueError:
                    continue
                options[tag] = value
                continue
            if c==header:
                header_line = line.rstrip().split(sep)
                c+=1
                continue
            c+=1
            line = line.rstrip().split(sep)
            construct_id = line[0]
            vals = {header_line[i]:line[i] for i in range(1,len(line))}
            if construct_id in ld:
                logging.error("{} {}".format(
                    "Duplicate construct ids found in library defintion:",
                    construct_id))
                raise RuntimeError()
            ld[construct_id] = vals
    return ld

# Tagging functions
def levenshtein_distance(observed, expected, allow_iupac=True):
    """ 
    Compute the levenshtein distance allowing for IUPAC encoding wildcards using
    the Wagner-Fischer algorithm.  
    
    Args:
        observed (str) - the observed barcode from sequencesin ACGT space
        expected (str) - the structure of the barcode in IUPAC format
        allow_iupac (bool) - True if IUPAC format is allowed for matching, False otherwise
    Returns:
        edit_distance (int) - The number of single character edits required to transform 
            the observed sequence into the expected sequence.

    Example:
        > levenshtein_distance("ACTG", "WSWS", allow_iupac=True)
        0
        > levenshtein_distance("ACTT", "WSWS", allow_iupac=True)
        1
        > levenshtein_distance("GACT", "WSWS", allow_iupac=True)
        2

    """
    if allow_iupac:
        translation = {"A":"A",
            "C":"C",
            "G":"G",
            "T":"T",
            "R":"AG",
            "Y":"CT",
            "S":"GC",
            "W":"AT",
            "K":"GT",
            "M":"AC",
            "B":"CGT",
            "D":"AGT",
            "H":"ACT",
            "V":"ACG",
            "N":"ACGT"} 
    else:
        translation = {"A":"A",
            "C":"C",
            "G":"G",
            "T":"T"}

    # get sizes and init an array to store edit distance counts
    lx = len(observed)+1
    ly = len(expected)+1
    m = np.zeros((lx,ly))
    # prepoulate edges
    m[:,0] = range(lx)
    m[0,:] = range(ly)
    for x in range(1,lx):
        for y in range(1,ly):
            if observed[x-1] in translation[expected[y-1]]:
                m[x,y] = m[x-1,y-1] # match, no penalty
            else:
                a = m[x-1,y]+1 # deletion
                b = m[x,y-1]+1 # insertion
                c = m[x-1,y-1]+1 # mismatch
                m[x,y] = min(a,b,c) # take the best one
    edit_distance = int(m[lx-1,ly-1]) 
    return edit_distance 

def hamming_distance(a, expected, allow_iupac=True):
    """ 
    Compute the hamming distance allowing for IUPAC encoding wildcards 
    
    Args:
        observed (str) - the observed barcode from sequencesin ACGT space
        expected (str) - the structure of the barcode in IUPAC format
        allow_iupac (bool) - True if IUPAC format is allowed for matching, False otherwise
    Returns:
        edit_distance (int) - The number number of mismatches 

    Example:
        > hamming_distance("ACTG", "WSWS", allow_iupac=True)
        0
        > hamming_distance("ACTT", "WSWS", allow_iupac=True)
        1
        > hamming_distance("GACT", "WSWS", allow_iupac=True)
        4

    """
    if allow_iupac:
        translation = {"A":"A",
            "C":"C",
            "G":"G",
            "T":"T",
            "R":"AG",
            "Y":"CT",
            "S":"GC",
            "W":"AT",
            "K":"GT",
            "M":"AC",
            "B":"CGT",
            "D":"AGT",
            "H":"ACT",
            "V":"ACG",
            "N":"ACGT"} 
    else:
        translation = {"A":"A",
            "C":"C",
            "G":"G",
            "T":"T"}


    c = 0
    for i in range(len(a)):
        try:
            if a[i] in translation[expected[i]]:
                continue
            else:
                c+=1
        except IndexError:
            c+=len(a)-i
    return c

def extract_subsequence(read, start = 30, end = 50, count_ed=False):
    """ Verify alignment against guide was good enough
        
        Args:
            read (pysam AlignedSegment) read to be analyzed
            guide_start (int): where the subsequence starts in the reference
            ends (int): where the subsequence ends in the reference
            extract (bool): True to return observed guide
        Return:
            guide_edit_distance (int) - number of one-nucleotide edits needed to transform the guide string into reference
            
    """
    # check if the read is rev_comp, if so the start and ends are calculated from the end of the rea
    ap = read.get_aligned_pairs(with_seq=True)
    if count_ed:
        matches = 0 
        mismatches = [0,0,0] # mismatches, insertions, deletions
    in_region=False
    i = 0
    subseq = ''
    while True:
        # get pair

        p = ap[i]
        i+=1     
        if p[1] is None:
            if in_region:
                if count_ed:
                    mismatches[1]+=1 # add insertion
                subseq+=read.query_sequence[p[0]]
            continue
        # only process the expected guide locations
        if p[1]<start:
            continue
        if p[1]==start:
            in_region=True
        if p[1]==end:
            break
        if in_region:
            # catalog mutations
            if p[0] is None:
                if count_ed:
                    mismatches[2]+=1 # add deletion to mismatches
                continue
            else:
                subseq+=read.query_sequence[p[0]]
                if count_ed:
                    if p[2].islower():
                        mismatches[0]+=1 # add mismatches
                    else:
                        matches += 1 # matches
        if i==len(ap):
            # never saw the end of the guide
            # so add the rest as deletions
            #logging.error(read.qname)
            if count_ed:
                mismatches[2] += end - start - sum(mismatches)  - matches
            break
    
    if subseq == '':
        subseq = None

    if count_ed:
        guide_edit_distance = sum(mismatches)
        return subseq, guide_edit_distance
    return subseq

def add_tags(read, guide_start=None, guide_length=None, expected_barcode=None,
    barcode_start=None, flag=None):
    """ Add tags to reads
        Count barcode edit distance

        Args:
            read (pysam AlignedSegment) - read to be analyzed
            guide_start (int): where the guide starts in the reference sequence
            guide_length (int): the guide length
            expected_barcode (str): the IPUAC structure of the barcode
            barcode_start (int): the position in the reference of the first base of the barcode
            flag (int): bitwise flag to add to the read (see sam format specification)
        Return:
            read (pysam AlignedSegment) - read with tags added            
    """
    if guide_start:
        #ged, seq = get_guide_edit_distance(read, start = guide_start, length=guide_length, extract=True)
        guide, ged = extract_subsequence(read, start = guide_start, end = guide_start+guide_length, count_ed=True)
    if expected_barcode:
        barcode = extract_subsequence(read, start = barcode_start, end = barcode_start + len(expected_barcode), count_ed=False)
        if barcode is not None:
            #bed = levenshtein_distance(barcode, expected_barcode)
            bed = hamming_distance(barcode, expected_barcode)

    # add tags to read
    # set guide edit distance
    if guide_start:
        read.set_tag('YG', guide, value_type="Z")
        read.set_tag('YH', ged, value_type="i")
    if expected_barcode and barcode:
        read.set_tag("YB", barcode, value_type = "Z")
        read.set_tag("YC", bed, value_type = "i")

    # add bitwise flags
    if flag:
        read.flag+=flag
    return read

# main functionality
def build_guide_reference(library=None,
    g_5p_r1=None, g_3p_r1=None, g_5p_r2=None, g_3p_r2=None,
    reference_fasta=None, tmp_dir = "./"):
    """
    Build expected guide sequences for aligner reference.

    Args: 
        library (dict): dictionary representation of the dual CRISPR guides 
            see: count.read_library_file
        g_5p_r1 (str): the 5' sequence of the construct upstream of the first
            guide.
        g_3p_r1 (str): the 3' sequence of the construct downstream of the first
            guide.
        g_5p_r2 (str): the 5' sequence of the construct upstream of the second
            guide.
        g_3p_r2 (str): the 3' sequence of the construct downstream of the second
            guide.
        reference_fasta (str): the name of the reference file
        tmp_dir (str): directory path to write the reference file to
    Returns:
        None
    Raises:
        RuntimeError if multiple sequences are observed for the same probe id
    """
    # check to see if at arguments are valid
    if library is None:
        raise RuntimeError('Must Supply a library defintion or file path.')
    if g_5p_r1 is None or g_3p_r1 is None:
        raise RuntimeError('Must specify read structrue for at least one read.')

    # read in library if a string is passed
    if isinstance(library,str):
        library = read_library_file(library)

    # build probe sequence dicts
    probe_a = defaultdict(list)
    probe_b = defaultdict(list)
    for construct in library:
        probe_a[library[construct]['probe_a_id']].append(
            library[construct]['probe_a_sequence'])

        try:
            probe_b[library[construct]['probe_b_id']].append(
                library[construct]['probe_b_sequence'])
        except KeyError:
            continue

    # write out 
    # check to see if an output reference file was passed
    if reference_fasta is None:
        # if none is passed then simply write to 
        # reference_fasta = os.path.join(tmp_dir,
        #                                 os.path.splitext(
        #                                     os.path.basename(library_path)
        #                                 )[0])
        reference_fasta = os.path.join(tmp_dir, 'reference.fa')

    with open(reference_fasta, 'w') as handle:

        for probe_id in sorted(probe_a):
            new_id = "{}_A".format(probe_id)
            seq = list(set(probe_a[probe_id]))
            if len(seq)>1:
                raise RuntimeError("{} {}".format(
                    "Probe sequences are not identical for",
                    probe
                    ))
            else:
                guide_sequence = seq[0]
            ref = g_5p_r1+guide_sequence+g_3p_r1
            handle.write(">{}\n".format(new_id))
            handle.write("{}\n".format(ref.upper()))             

        
        for probe_id in sorted(probe_b):
            new_id = "{}_B".format(probe_id)
            seq = list(set(probe_b[probe_id]))
            if len(seq)>1:
                raise RuntimeError("{} {}".format(
                    "Probe sequences are not identical for",
                    probe
                    ))
            else:
                guide_sequence = seq[0]
            ref = g_5p_r2+guide_sequence+g_3p_r2
            handle.write(">{}\n".format(new_id))
            handle.write("{}\n".format(ref.upper()))    
    return 0

def count_good_constructs(bam_file_path, 
    library = None,
    guide_edit_threshold = 2,
    sample = "Count",
    output_counts_path = "/dev/stdout"):
    """ Count construct pairs in aligned bam file.

    Args:
        bam_file_path (str): the bam file path containing the reads to count
        library (dict): dictionary representation of the dual CRISPR guides 
            see: count.read_library_file
        guide_edit_threshold (int): the maximum number of single base pair edits
            to allow in the guide region of aligned reads.
        output_counts_path (str): the path to write the construct counts to.
            Default is standard out.
    Return:
        None
    """
    # set up logger
    log = logging.getLogger()

    # check to make sure directory is writable
    if not os.path.exists(os.path.dirname(output_counts_path)):
        log.error("Directory of output path does not exists.")
        raise RuntimeError("Bad parameters.")
 

    if bam_file_path.endswith('bam'):
        method = "rb"
    else:
        method = "r"

    # open bam handle
    bam = pysam.AlignmentFile(bam_file_path, method)

    # get total read count and paired status
    log.info("Getting total read counts.")
    total_count, paired = get_fragment_count(bam_file_path)
    log.debug("Total Frags: {}. Is Paired: {}".format(total_count, paired))

    # initialize counters
    construct_counter = defaultdict(int)
    read_count = 0
    reads_considered = 0
    constructs_recognized = 0
    constructs_unrecognized = 0 
    guides_recognized = 0
    guides_unrecoginzed = 0
    valid_constructs = 0

    for read, mate in mate_pair_bam_reader(bam_file_path, paired=paired):

        read_count += 1
        # do some qc checks 
        if (read.is_unmapped) | (read.is_duplicate) | (paired and mate.is_unmapped) | (paired and not read.is_read1):
            continue

        # count the read
        reads_considered += 1 

        # get guide edit distance
        try:
            guide_ed_read1 = read.get_tag("YH")
        except KeyError:
            # no edit distance tag, assume bad
            continue

        # TODO: remove the positional tags
        rname = read.reference_name.rstrip("_A")
        #rname = read.reference_name

        # get construct
        if paired:
            mname = mate.reference_name.rstrip("_B")
            #mname = mate.reference_name
            construct = get_construct_divider().join([rname, mname])
            # get mate edit distance
            try:
                guide_ed_read2 = mate.get_tag("YH")
            except KeyError:
                # no edit distance tag, assume bad
                continue
        else:
            construct = rname
            guide_ed_read2 = None


        # Tabulate construct counts
        # count good guides
        # were the guides recognized according to allowed edit distance
        if (guide_ed_read1<=guide_edit_threshold):
            guides_recognized += 1
            if paired:
                if (guide_ed_read2<=guide_edit_threshold):
                    guides_recognized += 1
                    valid_constructs += 1
                    # report valid constructs
                    if library is not None:
                        if construct in library:
                            constructs_recognized+=1
                        else:
                            constructs_unrecognized+=1

                    # add to construct counter
                    construct_counter[construct]+=1
                else:
                    guides_unrecoginzed += 1

            else:
                construct_counter[construct]+=1
                valid_constructs += 1
                # report valid constructs
                if library is not None:
                    if construct in library:
                        constructs_recognized+=1
                    else:
                        constructs_unrecognized+=1
        else:
            guides_unrecoginzed += 1
            if paired:
                if (guide_ed_read2<=guide_edit_threshold):
                    guides_recognized += 1
                else:
                    guides_unrecoginzed += 1 


    if read_count==0:
        pct_valid = 0
        logging.warning("Dif ")
    else:
        pct_valid = valid_constructs/read_count * 100

    log.info("Found {0} passing constructs out of {1} reads. {2:.2f}%.".format(valid_constructs, read_count, pct_valid ))
    log.info("Writing outputs.")

    # finally write out counts and barcode paths
    with open(output_counts_path, 'w') as handle:
        handle.write("#Total Fragments: {}\n".format(total_count))
        handle.write("#Fragments Considered: {}\n".format(reads_considered))
        handle.write("#Valid Constructs: {}\n".format(valid_constructs))
        handle.write("#Guides Passing: {}\n".format(guides_recognized))
        handle.write("#Guides Failing: {}\n".format(guides_unrecoginzed))
        if library is not None:
            handle.write("#Constructs Recognized: {}\n".format(constructs_recognized))
            handle.write("#Constructs Unrecognized: {}\n".format(constructs_unrecognized))

        # write stats
        if library is not None:
            # use library to add write out extra columns
            # only write constructs we are interested in
            if paired: 
                header = ["construct_id",
                            "target_a_id","probe_a_id",
                            "target_b_id","probe_b_id",
                            sample]
            else:
                header = ["construct_id",
                            "target_a_id","probe_a_id",
                            sample]   
            out_line = "\t".join(header)
            handle.write("{}\n".format(out_line))
            for construct in sorted(library):
                info = library[construct]
                count = construct_counter[construct]
                to_write = [construct] + [info[i] for i in header[1:-1]] + [str(count)]
                out_line = "\t".join(to_write)
                try:
                    handle.write("{}\n".format(out_line))
                except UnicodeEncodeError:
                    print(out_line)
                    sys.exit()
        else:
            header = ["construct_id","count"]
            out_line = "\t".join(header)
            handle.write("{}\n".format(out_line))
            for construct in sorted(construct_counter):
                count = construct_counter[construct]
                to_write = [construct, str(count)]
                out_line = "\t".join(to_write)
                handle.write("{}\n".format(out_line))


    return 0

def count_good_constructs_and_barcodes(bam_file_path, 
    library = None,
    guide_edit_threshold = 2,
    barcode_edit_threshold = 0,
    sample = "Count",
    output_counts_path = "/dev/stdout",
    output_barcodes_path = None):
    """ Count construct pairs and valid barcodes.

    Args:
        bam_file_path (str): the bam file path containing the reads to count
        library (dict): dictionary representation of the dual CRISPR guides 
            see: count.read_library_file
        guide_edit_threshold (int): the maximum number of single base pair edits
            to allow in the guide region of aligned reads.
        barcode_edit_threshold (int): the maximum number of single base pair
            edits to allow in the barcode region of aligned reads.
        output_counts_path (str): the path to write the construct counts to.
            Default is standard out.
        output_barcodes_path (str): the path to write the observed barcode -
            construct mapping to.  Optional, default is to not write out
            barcodes.
    Return:
        None
    """

    # set up logger
    log = logging.getLogger()

    if bam_file_path.endswith('bam'):
        method = "rb"
    else:
        method = "r"

    # open bam handle
    bam = pysam.AlignmentFile(bam_file_path, method)

    # get total read count and paired status
    log.info("Getting total read counts.")
    total_count, paired = get_fragment_count(bam_file_path)
    log.debug("Total Frags: {}. Is Paired: {}".format(total_count, paired))

    # initialize counters
    construct_counter = defaultdict(int)
    observed_barcodes = defaultdict(lambda: defaultdict(int))
    read_count = 0
    reads_considered = 0
    constructs_recognized = 0
    constructs_unrecognized = 0 
    guides_recognized = 0
    guides_unrecoginzed = 0
    valid_barcodes = 0
    invalid_barcodes = 0
    barcode_assigned = 0
    valid_constructs = 0

    for read, mate in mate_pair_bam_reader(bam_file_path, paired=paired):

        read_count += 1
        # do some qc checks 
        if  (read.is_unmapped) | \
            (read.is_duplicate) | \
            (paired and mate.is_unmapped) | \
            (paired and not read.is_read1):
            continue

        # count the read
        reads_considered += 1 

        # get guide edit distance
        try:
            guide_ed_read1 = read.get_tag("YH")
        except KeyError:
            # no edit distance tag, assume bad
            continue


        rname = read.reference_name.rstrip("_A")
        #rname = read.reference_name

        # get construct
        if paired:
            mname = mate.reference_name.rstrip("_B")
            #mname = mate.reference_name
            construct = get_construct_divider().join([rname, mname])
            # get mate edit distance
            try:
                guide_ed_read2 = mate.get_tag("YH")
            except KeyError:
                # no edit distance tag, assume bad
                continue
        else:
            construct = rname
            guide_ed_read2 = None

        # TODO: fix this for no barcode reads
        # check barcode edit distance
        try:
            barcode = read.get_tag("YB")
            barcode_distance = read.get_tag("YC")
        except KeyError:
            # no barcode read, assume bad 
            log.debug("No barcode tag found in read: {}".format(read.qname))
            barcode = None 
            barcode_distance = barcode_edit_threshold+1 # one more than threshold so it never passes
            pass

        # Tabulate construct counts
        # count good guides
        # were the guides recognized according to allowed edit distance
        if (guide_ed_read1<=guide_edit_threshold):
            guides_recognized += 1
            if paired:
                if (guide_ed_read2<=guide_edit_threshold):
                    guides_recognized += 1
                    valid_constructs += 1
                    # report valid constructs
                    if library is not None:
                        if construct in library:
                            constructs_recognized+=1
                        else:
                            constructs_unrecognized+=1
                    # check the barcode
                    if (barcode_distance <= barcode_edit_threshold):
                        observed_barcodes[barcode][construct] += 1
                        barcode_assigned += 1
                        valid_barcodes += 1
                    else:
                        invalid_barcodes += 1

                    # add to construct counter
                    construct_counter[construct]+=1
                else:
                    guides_unrecoginzed += 1

            else:
                construct_counter[construct]+=1
                valid_constructs += 1
                # report valid constructs
                if library is not None:
                    if construct in library:
                        constructs_recognized+=1
                    else:
                        constructs_unrecognized+=1

                # check the barcode
                if (barcode_distance <= barcode_edit_threshold):
                    observed_barcodes[barcode][construct] += 1
                    barcode_assigned += 1
                    valid_barcodes += 1
                else:
                    invalid_barcodes += 1
        else:
            guides_unrecoginzed += 1
            if paired:
                if (guide_ed_read2<=guide_edit_threshold):
                    guides_recognized += 1
                else:
                    guides_unrecoginzed += 1 

            # check the barcode
            if (barcode_distance <= barcode_edit_threshold):
                observed_barcodes[barcode][construct] += 1
                #barcode_assigned += 1 # guides are bad so dont count
                valid_barcodes += 1
            else:
                invalid_barcodes += 1

    if read_count==0:
        pct_valid = 0
    else:
        pct_valid = valid_constructs/read_count * 100

    log.info("Found {0} passing constructs out of {1} reads. {2:.2f}%.".format(
        valid_constructs,
        read_count,
        valid_constructs/read_count*100))
    log.info("Writing outputs.")

    # finally write out counts and barcode paths
    with open(output_counts_path, 'w') as handle:
        handle.write("#Total Fragments: {}\n".format(total_count))
        handle.write("#Fragments Considered: {}\n".format(reads_considered))
        handle.write("#Guides Passing: {}\n".format(guides_recognized))
        handle.write("#Guides Failing: {}\n".format(guides_unrecoginzed))
        if library is not None:
            handle.write("#Constructs Recognized: {}\n".format(
                constructs_recognized))
            handle.write("#Constructs Unrecognized: {}\n".format(
                constructs_unrecognized))
        if valid_barcodes > 0:
            handle.write("#Barcodes Passing: {}\n".format(
                valid_barcodes))
            handle.write("#Barcodes Failing: {}\n".format(
                invalid_barcodes))
            handle.write("#Barcodes Assigned: {}\n".format(
                barcode_assigned))

        # write stats
        if library is not None:
            # use library to add write out extra columns
            # only write constructs we are interested in
            # TODO: this is hardcoded to two guides ... can we do singles?
            if paired: 
                header = ["construct_id",
                            "target_a_id","probe_a_id",
                            "target_b_id","probe_b_id",
                            sample]
            else:
                header = ["construct_id",
                            "target_a_id","probe_a_id",
                            sample]
            print(paired)     
            print(header)
            out_line = "\t".join(header)
            handle.write("{}\n".format(out_line))
            for construct in sorted(library):
                info = library[construct]
                count = construct_counter[construct]
                to_write = [construct] 
                to_write += [info[i] for i in header[1:-1]] 
                to_write += [str(count)]
                out_line = "\t".join(to_write)
                try:
                    handle.write("{}\n".format(out_line))
                except UnicodeEncodeError:
                    print(out_line)
                    sys.exit()
        else:
            header = ["construct_id","count"]
            out_line = "\t".join(header)
            handle.write("{}\n".format(out_line))
            for construct in sorted(construct_counter):
                count = construct_counter[construct]
                to_write = [construct, str(count)]
                out_line = "\t".join(to_write)
                handle.write("{}\n".format(out_line))

    if output_barcodes_path:
        # build ambigous barcodes path
        ambiguous = "{}.ambiguous.txt".format(
            os.path.splitext(output_barcodes_path)[0])
        with open(output_barcodes_path, 'w') as handle,\
             open(ambiguous, 'w') as amb:
            handle.write("barcode\tconstruct_id\tcount\n")
            amb.write("barcode\tconstruct_id\tcount\n") 
            for barcode, vals in observed_barcodes.items():
                constructs = [(i[0], i[1]) for i in vals.items()]
                if len(constructs)>1: 
                    log.debug("Abgious Barcode: {}. Maps to {} constructs".\
                        format(barcode, len(constructs)))
                    output = "{}\t".format(barcode) 
                    for construct in constructs:
                        output += "{}\t{}\t".format(construct[0],construct[1])
                    amb.write("{}\n".format(output))
                else:
                    # split out information
                    construct = constructs[0][0] # name of constructs
                    count = constructs[0][1]
                    to_write = [barcode, construct, str(count)]
                    out_line = "\t".join(to_write)
                    handle.write("{}\n".format(out_line))

    return 0

