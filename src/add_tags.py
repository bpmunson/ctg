#!/usr/bin/env python
import sys
import pysam
import argparse
import regex as re

# Tagging functions
def levenshtein_regex(a, expected):
    """
    Calculate the edit distance of a sequence from the expected sequence
    structure.
    Args:
        a (str): the sequence to test
        expected (str): the expected fromat of a in IUPAC format.
    Returns:
        edit_distance (int): the minimum number of 1 base pair changes required
        to make the sequence into the expected sequence structure.
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
    edit_distance = sum(h.fuzzy_counts)
    return edit_distance

def get_guide_edit_distance(read, start = 30, length = 20, extract=False):
    """ Calculate the edit distance in a similar vein to the sam format
        'NM tag', for a subregion on the reference corresponding to the guide.

        Args:
            read (pysam AlignedSegment): the read segment to inspect
            start (int): where the first base of the subregion in the reference
            length (int): where the length of the subregion
            extract (bool): return the query subsequence as well as the 
                edit distance to the reference
        Return:
            guide_edit_distance (int): number of one-nucleotide edits
                needed to transform the subsequence into reference
            subseq (str): the query subsequence (Optional)
        Raises:
            IndexError: If the subsequence extends beyond the reference.
    """
    # check if the read is rev_comp, if so the start and
    # ends are calculated from the end of the read
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
            logging.error("Index out of range: {} {}".format(p, i))
            raise IndexError()
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

def extract_barcode(read, start=175, length=30):
    """ Extract the UMI barcode (randomer) from a specified position in 
        the reference

        Args:
            read (pysam AlignedSegment): the read segment to inspect
            start (int): where the first base of the subregion in the reference
            length (int): where the length of the subregion
        Return:
            barcode (str): the query sequence in the read, corresponding to 
                the subregion of the reference sequence specifed.
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
        if ap[i][1] == start:
            ap_start = i
        elif ap[i][1] == start+length:
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

def add_tags(read,
            guide_start=None, guide_length=None,
            expected_barcode=None, barcode_start=None,
            flag=None):
    """ Add sam/bam tags to a pysam read alignment. 
        Optionally add information about the edit distance in the guide region,
        the randomer barcode, and mate information stored in the bitwise flag.

        Args:
            read (pysam AlignedSegment): the read segment to inspect
            guide_start (int): where the first base of the subregion in the
                reference where the guide occurs
            guide_length (int): where the length of guide in the reference
            expected_barcode (str): the expected format of the barcode in
                IUPAC format.
            barcode_start (int): the start position of the barcode in 
                the reference
            flag (int): bitwise flag to add to the the read
        Return:
            read (pysam AlignedSegment): the read segment with modififications
    """
    if guide_start:
        ged, seq = get_guide_edit_distance( read,
                                            start = guide_start,
                                            length=guide_length,
                                            extract=True)
    if expected_barcode:
        barcode = extract_barcode(read,
                                start = barcode_start,
                                length = len(expected_barcode))
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

def main():
    """ Main - Command line functionality to add reads to streamed sam file.

    Input is a piped sam file.
    Output is a piped sam file with tags added to the reads.
    """
    parser = argparse.ArgumentParser(description="Add tags to aligned reads.")
    parser.add_argument("--flag",
        action="store", default=None, type=int,
        help="Bitwise flags to add to all reads.")
    parser.add_argument("--expected_barcode",
        action="store", type=str, default=None, required=False,
        help="Position of the guide in the read.")
    parser.add_argument("--barcode_start",
        action="store",type=int,default=None,required=False,
        help="Position of barcode in read")
    parser.add_argument("--guide_start",
        action="store", type=int,
        help="Position of the guide in the read.")
    parser.add_argument("--guide_length",
        action="store",type=int,
        help="Position of the guide in the read.")
    args = parser.parse_args()

    bam = pysam.AlignmentFile("-", "r")

    # get header
    header= bam.header.to_dict()
    # add new process group tag to header
    new_pg = {'CL':' '.join(sys.argv), 'ID':'CTG', 'PN':'add_tags', 'VN':0.1}
    header['PG'].append(new_pg)
    # open output bam handle
    output = pysam.AlignmentFile("-", "w", header=header)

    # loop through reads 
    while True:
        try:
            read = next(bam)
        except StopIteration:
            break

        # if a read is unmapped just add the tags 
        if read.is_unmapped:
            read = add_tags(read,
                            flag=args.flag)
        else:
            read = add_tags(read,
                            args.guide_start,
                            args.guide_length
                            args.expected_barcode,
                            args.barcode_start,
                            args.flag)
        try:
            output.write(read)
        except OSError:
            # this is to handle keyboard interupts
            # or pipeling to premature stop (eg: | head)
            break

    # close the handles
    bam.close()
    output.close()

if __name__ == "__main__":
    main()