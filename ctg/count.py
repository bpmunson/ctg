
"""
Docstring
"""

import os
import pysam
import logging
import regex as re
from collections import defaultdict

import pandas as pd

_construct_divider = "__"
_probe_divider = "-"

# get complement of a sequence
comp = lambda x: x.translate(str.maketrans('ACGT','TGCA'))

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
        RunTimeError if duplicate constructs were seen in library 

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
                raise RuntimeError("{} {}".format(
                    "Duplicate construct ids found in library defintion:",
                    construct_id))
            ld[construct_id] = vals
    return ld

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

# Tagging functions
def levenshtein_distance(observed, expected, allow_iupac=True):
    """ 
    Compute the levenshtein distance allowing for IUPAC encoding wildcards
    Args:
        observed (str) - the observed barcode from sequencesin ACGT space

    # seq is the observed barcode sequence
    # expected is the IUPAC enocded structure
    """
    if allow_iupac:
        IUPAC = {"A":"A",
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
        IUPAC = {"A":"A",
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
            a = m[x-1,y]+1 # insertion
            b = m[x,y-1]+1 # deletion
            if (observed[x-1] in IUPAC[expected[y-1]]):
                c = m[x-1,y-1] # match
            else:
                c = m[x-1,y-1]+1 # mismatch
            m[x,y] = min(a,b,c) # take the best one
    return (m[lx-1,ly-1]) 

def levenshtein_regex(a, expected):
    """ Return the fuzzy counts of a regex match to the expected barcode structure.

        Args:
            a (str): string to test
            expected (str) : [AC][GT][AC][GT] ....
        Return:
            edit distance
    """
    IUPAC = {"A":"A",
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
    # consturct regex string
    expected = "".join(["[{}]".format(IUPAC[b]) for b in expected])
    r = re.compile("({}){{e}}".format(expected))
    h = r.match(a)
    return sum(h.fuzzy_counts)

def get_guide_edit_distance(read, start = 30, length = 20, extract=False):
    """ Verify alignment against guide was good enough
        
        Args:
            read (pysam AlignedSegment) read to be analyzed
            guide_start (int): where the guide starts in the reference sequence
            length (int): where the guide ens in the reference sequence
            extract (bool): True to return observed guide sequence
        Return:
            guide_edit_distance (int) - number of one-nucleotide edits needed to transform the guide string into reference
            
    """
    # check if the read is rev_comp, if so the start and ends are calculated from the end of the rea
    ap = read.get_aligned_pairs(with_seq=True)
    matches = 0 
    mismatches = [0,0,0] # mismatches, insertions, deletions
    in_region=False
    i = 0
    if extract:
        subseq = ''
    while True:
        # get pair
        try:
            p = ap[i]
        except IndexError:
            print(i)
            print(ap)
            print(read.qname)
            sys.exit()
        i+=1     
        if p[1] is None:
            if in_region:
                mismatches[1]+=1 # add insertion
                if extract:
                    subseq+=read.query_sequence[p[0]]
            continue
        # only process the expected guide locations
        if p[1]<start:
            continue
        if p[1]==start:
            in_region=True
        if p[1]==start+length:
            break
        if in_region:
            # catalog mutations
            if p[0] is None:
                mismatches[2]+=1 # add deletion to mismatches
            # elif p[1] is None:
            #     mismatches[1]+=1 # add insertion
            elif p[2].islower():
                mismatches[0]+=1 # add mismatches
                if extract:
                    subseq+=read.query_sequence[p[0]]

            else:
                matches += 1 # matches
                if extract:
                    subseq+=read.query_sequence[p[0]]
        if i==len(ap):
            # never saw the end of the guide
            # so add the rest as deletions
            #logging.error(read.qname)
            
            mismatches[2]+= length - sum(mismatches)  - matches
            break
    
    guide_edit_distance = sum(mismatches)
    # return pass
    if extract:
        return guide_edit_distance, subseq 
    return guide_edit_distance

def extract_barcode(read, s=175, e=175+30):
    """ Extract random-mer barcode from read
        Todo: we shouldn't have to loop through read twice
        
        Args:
            pysam AlignedSegment
            barcode_start - where the guide starts in the reference sequence
            barcode_end - where the guide ens in the reference sequence
            threshold - the maximum hamming distance of the barcode to an expected one
        Return:
            barcode sequence            
    """
    #read_length = read.infer_query_length()
    ap = read.get_aligned_pairs(with_seq=True)
    barcode = ''
    # init
    ap_start, ap_end = None, None
    # find start in pairs
    for i in range(len(ap)):
        if ap[i][1] is None:
            continue
        if ap[i][1] == s:
            ap_start = i
        elif ap[i][1] == e:
            ap_end = i

    if ap_start is None:
        # never saw start in reference
        return None
    if ap_end is None:
        # never saw the end in the refernce
        return None

    i = ap_start
    while True:
        if i==ap_end:
            break
        if ap[i][0]==None: # deletion
            ap_end+=1
            i+=1
            continue
        if ap[i][0]==None: # insertion
            ap_end-=1
        barcode += read.query_sequence[ap[i][0]]
        i+=1
        
    if barcode == "":
        return None
    return barcode 

def add_tags(read, guide_start=None, guide_length=None, expected_barcode=None, barcode_start=None, flag=None):
    """ Add tags to reads
        Count barcode edit distance
    """
    if guide_start:
        ged, seq = get_guide_edit_distance(read, start = guide_start, length=guide_length, extract=True)
    if expected_barcode:
        barcode = extract_barcode(read, s = barcode_start, e = barcode_start + len(expected_barcode))
        if barcode is not None:
            bed = levenshtein_regex(barcode, expected_barcode)

    # add tags to read
    # set guide edit distance
    if guide_start:
        read.set_tag('YG', seq, value_type="Z")
        read.set_tag('YH', ged, value_type="i")
    if expected_barcode and barcode:
        read.set_tag("YB", barcode, value_type = "Z")
        read.set_tag("YC", bed, value_type = "i")

    # add bitwise flags
    if flag:
        read.flag+=flag
    return read

# main functionality
def build_guide_reference(library,
    g_5p_r1, g_3p_r1, g_5p_r2, g_3p_r2,
    reference_fasta=None, tmp_dir = "./"):
    """
    Build expected guide sequences for aligner reference

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
        RunTimeError if multiple sequences are observed for the same probe id
    """

    # read in library if a string is passed
    if isinstance(library,str):
        library, options = read_library_file(library_path)

    # build probe sequence dicts
    probe_a = defaultdict(list)
    probe_b = defaultdict(list)
    for construct in library:
        probe_a[library[construct]['probe_a_id']].append(
            library[construct]['probe_a_sequence'])
        probe_b[library[construct]['probe_b_id']].append(
            library[construct]['probe_b_sequence'])

    # write out 
    # check to see if an output reference file was passed
    if reference_fasta is None:
        # if none is passed then simply write to 
        reference_fasta = os.path.join(tmp_dir,
                                        os.path.splitext(
                                            os.path.basename(library_path)
                                        )[0])
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
    output_counts_path = "/dev/stdout"):
    """ Count construct pairs and assign barcodes
    """

    # set up logger
    log = logging.getLogger()

    if bam_file_path.endswith('bam'):
        method = "rb"
    else:
        method = "r"

    # open bam handle
    bam = pysam.AlignmentFile(bam_file_path)

    # get total read count and paired status
    log.info("Getting total read counts.")
    total_count, paired = get_fragment_count(bam_file_path)
    log.debug("Total Frags: {}. Paired: {}".format(total_count, paired))

    # initialize counters
    construct_counter = defaultdict(int)
    observed_barcodes = defaultdict(lambda: defaultdict(int))
    read_count = 0
    reads_considered = 0
    constructs_recognized = 0
    constructs_unrecognized = 0 
    guides_recognized = 0
    guides_unrecoginzed = 0
    valid_constructs = 0

    for read, mate in mate_pair_bam_reader(bam_file_path):

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
            if (guide_ed_read2<=guide_edit_threshold):
                guides_recognized += 1
            else:
                guides_unrecoginzed += 1 



    log.info("Found {0} passing constructs out of {1} reads. {2:.2f}%.".format(valid_constructs, read_count, valid_constructs/read_count*100 ))
    log.info("Writing outputs.")

    # finally write out counts and barcode paths
    with open(output_counts_path, 'w') as handle:
        handle.write("#Total Fragments: {}\n".format(total_count))
        handle.write("#Fragments Considered: {}\n".format(reads_considered))
        handle.write("#Guides Passing: {}\n".format(guides_recognized))
        handle.write("#Guides Failing: {}\n".format(guides_unrecoginzed))
        if library is not None:
            handle.write("#Constructs Recognized: {}\n".format(constructs_recognized))
            handle.write("#Constructs Unrecognized: {}\n".format(constructs_unrecognized))

        # write stats
        if library is not None:
            # use library to add write out extra columns
            # only write constructs we are interested in
            header = ["construct_id","target_a_id","probe_a_id","target_b_id","probe_b_id","count"]
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
    bam = pysam.AlignmentFile(bam_file_path)

    # get total read count and paired status
    log.info("Getting total read counts.")
    total_count, paired = get_fragment_count(bam_file_path)
    log.debug("Total Frags: {}. Paired: {}".format(total_count, paired))

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

    for read, mate in mate_pair_bam_reader(bam_file_path):

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
            header = ["construct_id",
                        "target_a_id","probe_a_id",
                        "target_b_id","probe_b_id",
                        "count"]
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

def aggregate_counts(libary, counts_files,
    output_file = "/dev/null", names = None):
    """ Combine counts files by constrcut
    """
    counts_files = ['counts_1.txt', 'counts_2.txt']
    library = "library_defintion/library_batch_1_9.csv"
    names = [1,2]
    #combined_counts_dictionary = defaultdict(lambda: defaultdict(int))
    # use pandas to merge together
    dfs = [pd.read_csv(file, sep="\t", header=0, index_col=0) for file in counts_files]
    cc = pd.concat(dfs, axis=1)
    cc.index.name="construct_id"
    cc.columns = names

    # get library
    ld = pd.read_csv(library, sep="\t", header=0)

    # subset library definition
    info = ld[['target_a_id','probe_a_id','target_b_id','probe_b_id']]
    # build pair names
    construct_id = list(ld[['probe_a_id','probe_b_id']].apply(lambda x: '__'.join(x), axis=1))
    probe_pair_id = list(ld[['probe_a_id','probe_b_id']].apply(lambda x: '__'.join(sorted(x)),axis=1))
    idx = list(info[['probe_a_id','probe_b_id']].apply(lambda x: "{}_A__{}_B".format(x[0],x[1]),axis=1))
    # add info to data frame
    info.loc[:,'construct_id'] = construct_id
    info.loc[:,'probe_pair_id'] = probe_pair_id
    # set index 
    info.loc[:,'idx'] = idx
    info.set_index('idx',inplace=True)
    # reorder 
    info = info[['construct_id','target_a_id','probe_a_id','target_b_id','probe_b_id','probe_pair_id']]

    combined_counts = pd.merge(info, cc, left_index=True, right_index=True)
    combined_counts.to_csv(output_file, header=True, index=False, sep="\t")

    return 0

