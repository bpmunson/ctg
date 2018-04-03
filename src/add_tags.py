#!/usr/bin/env python
import sys
import pysam
import argparse
import regex as re
#from align import levenshtein_regex, extract_barcode, get_guide_edit_distance



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
            a: string to test
            expected : [AC][GT][AC][GT] ....
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
            pysam AlignedSegment
            guide_start - where the guide starts in the reference sequence
            guide_end - where the guide ens in the reference sequence
            threshold - the maximum mutations (insertions, deletions, substitions) allowed
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


def main():
    parser = argparse.ArgumentParser(description="Add tags to aligned reads.")
    parser.add_argument("--flag", action="store", default=None, type=int, help="Bitwise flags to add to all reads.")
    parser.add_argument("--expected_barcode", action="store", type=str, default=None,required=False, help="Position of the guide in the read.")
    parser.add_argument("--barcode_start",action="store",type=int,default=None,required=False,help="Position of barcode in read")
    parser.add_argument("--guide_start", action="store", type=int,help="Position of the guide in the read.")
    parser.add_argument("--guide_length", action="store",type=int, help="Position of the guide in the read.")
    #parser.add_argument("input_bam",action="store",help="input bam")
    #parser.add_argument("output_bam",action="store",help="output bam")
    args = parser.parse_args()

    bam = pysam.AlignmentFile("-", "r")

    # get header
    header= bam.header.to_dict()
    new_pg = {'CL':' '.join(sys.argv), 'ID':'CTG', 'PN':'add_tags', 'VN':0.1}
    header['PG'].append(new_pg)
    output = pysam.AlignmentFile("-", "w", header=header)

    while True:
        try:
            read = next(bam)
        except StopIteration:
            break

        if not read.is_unmapped:
            read = add_tags(read, args.guide_start, args.guide_length, args.expected_barcode, args.barcode_start, args.flag)
        else:
            read = add_tags(read, flag=args.flag)

        try:
            output.write(read)
        except OSError:
            break

    bam.close()
    output.close()


if __name__ == "__main__":
    main()