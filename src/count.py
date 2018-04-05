"""
General Pipeline

From library definition file and expected read structure, make bowtie2 reference
Align reads independently to guide reference
Merge Reads together and add mate information
Count good constructs while extracting barcodes  


guide_5p_r1: the 5' side of the first guide 
guide_3p_r1: the 3' side of the first guide 
guide_5p_r2: the 5' side of the second guide in the same orientation as the first)
guide_3p_r2: the 3' side of the second guide in the same orientation as the first)

ie if the construct is 
R1|--------------->|
5'-AAAA<Guide1>CCCCNNNNGGGG<Guide2>TTTT-3'
                      |<---------------|R2
guide_5p_r1=AAAA
guide_3p_r1=CCCC
guide_5p_r2=GGGG
guide_3p_r2=TTTT
"""

import os
import shutil
import time
import pysam
import logging
import tqdm

import pandas as pd


from collections import defaultdict

import align



_construct_divider = "__"
_probe_divider = "-"

def get_construct_divider():
    return _construct_divider

def get_probe_divider():
    return _probe_divider

# get complement of a sequence
comp = lambda x: x.translate(str.maketrans('ACGT','TGCA'))

# get reverse complement of seqeunce
rev_comp = lambda x: comp(x)[::-1]


class TqdmLoggingHandler(logging.Handler):
    """
    Taken from: https://stackoverflow.com/questions/38543506/change-logging-print-function-to-tqdm-write-so-logging-doesnt-interfere-wit/38739634#38739634

    ex: 

    import logging
    import time
    import tqdm
    import io


    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    logger.addHandler (TqdmLoggingHandler (level=logging.INFO))

    for x in tqdm.tqdm(range(10)):
        if x == 5:
            logging.info("Half way")
        time.sleep(.1)


    """
    def __init__ (self, level = logging.INFO):
        super (self.__class__, self).__init__(level)

    def emit (self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)      

class Counter():

    def __init__(self, library, fastq1,
        fastq2=None,
        input_bam=None,
        barcode=None, 
        barcode_location = 175,
        barcode_read = '1',
        guide_edit_threshold=2,
        barcode_edit_threshold =0,
        guide_5p_r1 = "TATATATCTTGTGGAAAGGACGAAACACCG",
        guide_3p_r1 = "GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGGCCACCTTAACACGCGATGATATTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCTATTACGAGCGCTTGGAT",
        guide_5p_r2 = "CTTGGAGAAAAGCCTTGTTTG",
        guide_3p_r2 = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG",
        guide_length_r1 = 20,
        guide_length_r2 = 20,
        output_counts = "/dev/stdout",
        output_bam = None,
        output_barcodes = None,
        sample = None,
        time_point = None,
        replicate = None,
        get_name_from_fastq = False,
        threads = 1):

        # read in library 
        self.library = read_library_file(library)
        
        # store args
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.barcode = barcode
        self.barcode_location = barcode_location
        self.barcode_read = barcode_read
        self.guide_5p_r1 = guide_5p_r1
        self.guide_3p_r1 = guide_3p_r1
        self.guide_5p_r2 = guide_5p_r2
        self.guide_3p_r2 = guide_3p_r2
        self.guide_length_r1 = guide_length_r1
        self.guide_length_r2 = guide_length_r2
        self.guide_edit_threshold = guide_edit_threshold
        self.barcode_edit_threshold = barcode_edit_threshold
        
        self.output_bam = output_bam
        self.output_counts = output_counts
        self.output_barcodes = output_barcodes

        self.threads = threads
        # get directory
        self.dir = os.getcwd()

        # make a temporary working directory 
        self.tmp_dir = os.path.join(self.dir,"ctg_tmp_{}".format(int(time.time())))
        os.makedirs(self.tmp_dir)

        self.output_summary_path = os.path.join(self.tmp_dir, "tmp.summary.csv")

        # store aligned_bam if supplied
        self.aligned_bam = input_bam


        # check the supplied parameters
        self._check_args()

    def _check_args(self):
        """ Verify arguments pass some logic checks
        """
        if self.barcode:
            if (self.barcode_read == "2"):
                if (self.fastq2 is None):
                    logging.error("Barcode was specified in read 2 but no read 2 fastq was passed.")
                    raise RuntimeError("Bad arguments")
                if self.barcode_location +len(self.barcode)> len(self.guide_5p_r2)+len(self.guide_3p_r2)+self.guide_length_r2:
                    logging.error("Barcode position extends beyond expected read length.")
                    raise RuntimeError("Bad arguments")
            else:
                if self.barcode_location +len(self.barcode)> len(self.guide_5p_r1)+len(self.guide_3p_r1)+self.guide_length_r1:
                    logging.error("Barcode position extends beyond expected read length.")
                    raise RuntimeError("Bad arguments")

    def _parse_fastq_name(self):
        """ Parse fastq name to get sample, time point, and replicate information
            Assume fastq is formatted as:
                <SAMPLE>_[d|t]<Timepoint>_<Replicate>_...fastq
        """
        name = os.path.splitext(os.path.basename(self.fastq1))[0]
        split = name.split("_")


        if len(split)<3:
            logging.error("Incorrect fastq name structure, reverting to default. ex: SAMPLE_d3_1_r1.fastq.")
            return None

        sample = split[0]
       
        timepoint = split[1]
        if timepoint[0].upper() not in ['T','D']:
            logging.error("Parsing fastq file name. Timepoint style string not in second field, reverting to default.  Please use correct file name. ex: SAMPLE_d3_1_r1.fastq.")
            return None
        try:
            replicate = split[2]
            _ = int(replicate)
        except ValueError:
            logging.warning("Parsing fastq file name. Third field doesn't appear to be a number. Using default.")
            #return None


        fid = "_".join(split[:3])
        return fid

    def build_reference(self):
        """ Build reference from library definition
        """
        if self.aligned_bam:
            logging.warning("An aligned bam was provided as an input, negating the need to build a reference.")
            return 0

        self.reference = os.path.join(self.tmp_dir, "guide_reference.fa")

        # parse library definition and write out guide refernece
        build_guide_reference(  self.library,
                                self.guide_5p_r1,
                                self.guide_3p_r1,
                                self.guide_5p_r2,
                                self.guide_3p_r2,
                                reference_fasta=self.reference,
                                tmp_dir=self.tmp_dir)

        # run bowtie build on the reference
        align.bowtie2_build(self.reference, self.reference, "--quiet")

        return 0

    def align_reads(self):
        """ Run bowtie to align to guide reference 
        """
        if self.aligned_bam:
            logging.warning("An aligned bam was provided as an input, negating the need to align.")
            return 0

        # align read 2
        if self.fastq2:
            r1_flags_to_add = 65
            r2_flags_to_add = 129
        else:
            r1_flags_to_add = 0
            r2_flags_to_add = None

        if self.barcode_read == "1":
            barcode=self.barcode
            barcode_location = self.barcode_location
        else:
            barcode = None
            barcode_location = None

        # align read 1
        bam1 = os.path.join(self.tmp_dir, "tmp.r1.bam")
        align.bowtie2_align(self.fastq1, self.reference, "--quiet", "--very-sensitive",
            bam=bam1,
            binary_flags_to_add=r1_flags_to_add,
            guide_start=len(self.guide_5p_r1), 
            guide_length=self.guide_length_r1,
            expected_barcode=barcode,
            barcode_start=barcode_location,
            threads=self.threads )


        if self.fastq2:
            if self.barcode_read == "2":
                barcode=self.barcode
                barcode_location = self.barcode_location
            else:
                barcode = None
                barcode_location = None
            # align read 2 
            bam2 = os.path.join(self.tmp_dir, "tmp.r2.bam")
            align.bowtie2_align(self.fastq2, self.reference, "--quiet", "--very-sensitive",
                bam=bam2,
                binary_flags_to_add=r2_flags_to_add,
                guide_start=len(self.guide_5p_r2), 
                guide_length=self.guide_length_r2,
                expected_barcode=barcode,
                barcode_start=barcode_location,
                threads=self.threads)

            # merge
            self.aligned_bam = os.path.join(self.tmp_dir, "tmp.bam")
            align.merge_fixmate(bam1, bam2, self.aligned_bam, tmp=os.path.join(self.tmp_dir,"sort"))

            # remove temporary read level bams
            os.remove(bam1)
            os.remove(bam2)

        else:
            # if no read 2, just keep read 1 aligned bam and store path
            self.aligned_bam = bam1

    def count(self):
        """ Run counting
        """
        log = logging.getLogger()
        if self.aligned_bam is None:
            raise RuntimeError("Must align reads before counting.")

        log.debug("Aligned bam: {}".format(self.aligned_bam))
        
        if (not os.path.exists(self.aligned_bam)):
            raise RuntimeError("Bam does not exists. Must first align fastqs.")

        log.info("Couting Aligned Constructs.")
        count_good_constructs(self.aligned_bam,
            library = self.library, 
            guide_edit_threshold = self.guide_edit_threshold,
            barcode_edit_threshold = self.barcode_edit_threshold,
            output_counts_path = self.output_counts,
            output_barcodes_path = self.output_barcodes)


        return

    def cleanup(self):
        """ Remove tmp directory """
        logging.info("Cleaning up.")
        # store the aligned bam if desired
        if self.output_bam is not None:
            shutil.move(self.aligned_bam, self.output_bam)

        # delete temporary tree
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

        try:
            shutl.rmdir(self.tmp_dir)
        except:
            pass

### Helpers


def read_1_callback(read):
    """ Function to determine if read is the first in a pair
        Args:
            pysam AlignedSegment
        Return:
            True/False
    """
    if read.is_read1:
        return True
    else:
        return False

def read_library_file(library, sep="\t", header=0, comment="#"):
    """ Open read in a library file
    format is:
    construct_id    target_a_id probe_a_id  probe_a_sequence    target_b_id probe_b_id  probe_b_sequence

    constructs_id is 
    """
    ld = dict()
    with open(library, 'r') as handle:
        c=0
        for line in handle:
            if line.startswith(comment):
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
                raise RuntimeError("Duplicate construct ids found in library defintion: {}".format(construct_id))
            ld[construct_id] = vals
    return ld

def mate_pair_bam_reader(bam_file_path, paired=True):
    """
    Read bam file mate pairs
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
    Get the number of fragments an aligned bam, checking for paired status first

    Returns the readcount and boolean indicating paired status
    """
    # check paired
    read_count = int(pysam.view("-c", "-f65", bam).rstrip())
    if read_count == 0:
        # get singles
        read_count = int(pysam.view("-c",bam).rstrip())
        return read_count, False
    else:
        return read_count, True

def build_guide_reference(library, g_5p_r1, g_3p_r1, g_5p_r2, g_3p_r2, reference_fasta=None, tmp_dir = "./"):
    """
    Build expected guide sequences for aligner reference

    Args: output reference
    """

    # read in library if a string is passed
    if isinstance(library,str):
        library = read_library_file(library_path)

    # build probe sequence dicts
    probe_a = defaultdict(list)
    probe_b = defaultdict(list)
    for construct in library:
        probe_a[library[construct]['probe_a_id']].append(library[construct]['probe_a_sequence'])
        probe_b[library[construct]['probe_b_id']].append(library[construct]['probe_b_sequence'])

    # write out 
    # check to see if an output reference file was passed
    if reference_fasta is None:
        # if none is passed then simply write to 
        reference_fasta = os.path.join(tmp_dir, os.path.splitext(os.path.basename(library_path))[0])
    with open(reference_fasta, 'w') as handle:

        for probe_id in sorted(probe_a):
            new_id = "{}_A".format(probe_id)
            seq = list(set(probe_a[probe_id]))
            if len(seq)>1:
                raise RuntimeError("Probe sequences are not identical for {}".format(probe))
            else:
                guide_sequence = seq[0]
            ref = g_5p_r1+guide_sequence+g_3p_r1
            handle.write(">{}\n".format(new_id))
            handle.write("{}\n".format(ref.upper()))             

        for probe_id in sorted(probe_b):
            new_id = "{}_B".format(probe_id)
            seq = list(set(probe_b[probe_id]))
            if len(seq)>1:
                raise RuntimeError("Probe sequences are not identical for {}".format(probe))
            else:
                guide_sequence = seq[0]
            ref = g_5p_r2+guide_sequence+g_3p_r2
            handle.write(">{}\n".format(new_id))
            handle.write("{}\n".format(ref.upper()))    
    return 0

def count_good_constructs(bam_file_path, 
    library = None,
    guide_edit_threshold = 2,
    minimum_alignment_score = -55.5,
    barcode_edit_threshold = 0,
    output_counts_path = "/dev/stdout",
    output_barcodes_path = None,
    show_progress_bar=True):
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

    if show_progress_bar == True:
        # get old handlers
        old_handlers = log.handlers
        # remove them
        log.handlers = []
        # add new handler
        log.addHandler(TqdmLoggingHandler(level=logging.INFO))

        pbar = tqdm.tqdm(total=total_count)
        m = 100

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
        if show_progress_bar:
            if read_count % m == 0:
                pbar.update(m) 
        
        # do some qc checks 
        if (read.is_unmapped) | (read.is_duplicate) | (paired and mate.is_unmapped) | (paired and not read.is_read1):
            continue

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
            barcode_ditance = barcode_edit_threshold+1 # one more than threshold so it never passes
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


    if show_progress_bar:
        pbar.close()
        # remove extra handler
        log.handlers = old_handlers

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
        if valid_barcodes > 0:
            handle.write("#Barcodes Passing: {}\n".format(valid_barcodes))
            handle.write("#Barcodes Failing: {}\n".format(invalid_barcodes))
            handle.write("#Barcodes Assigned: {}\n".format(barcode_assigned))

        # write stats
        if library is not None:
            # use library to add write out extra columns
            # only write constructs we are interested in
            header = ["construct_id","count","target_a_id","probe_a_id","target_b_id","probe_b_id"]
            out_line = "\t".join(header)
            handle.write("{}\n".format(out_line))
            for construct in sorted(library):
                info = library[construct]
                count = construct_counter[construct]
                to_write = [construct, str(count)] + [info[i] for i in header[2:]]
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
        ambiguous = "{}.ambiguous.txt".format(os.path.splitext(output_barcodes_path)[0])
        with open(output_barcodes_path, 'w') as handle, open(ambiguous, 'w') as amb:
            handle.write("barcode\tconstruct_id\tcount\n")
            amb.write("barcode\tconstruct_id\tcount\n") 
            for barcode, vals in observed_barcodes.items():
                constructs = [(i[0], i[1]) for i in vals.items()]
                if len(constructs)>1: 
                    log.debug("Abgious Barcode: {}. Maps to {} constructs".format(barcode, len(constructs)))
                    output = "{}\t".format(barcode) 
                    for construct in constructs:
                        output += "{}\t{}\t".format(construct[0],construct[1])
                    amb.write("{}\n".format(output))
                else:
                    # split out information
                    construct = constructs[0][0] # name of constructs
                    count = constructs[0][1]
                    to_write = [construct, str(count)]
                    out_line = "\t".join(to_write)
                    handle.write("{}\n".format(out_line))

    return 0

def aggregate_counts(libary, counts_files, output_file = "/dev/null", names = None):
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

