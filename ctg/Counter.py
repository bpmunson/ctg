import os
import time
import shutil
import logging
from ctg.core import count
from ctg.core import align

class Counter():
    """ CTG Counter: count dual knockout contructs from sequencing data.
    
    Args:
        library (str): Library defintion file in csv format.
        fastq1 (str): File with read 1 mates. Can be gzipped.
        guide_5p_r1 (str): Expected 5' end of read 1.
        guide_3p_r1 (str): Expected 3' end of read 1.
        guide_5p_r2 (str): Expected 5' end of read 2.
        guide_3p_r2 (str): Expected 3' end of read 2.
        threads (int): Number of threads to use.
        sample (str): Sample name to add to counts file.
        fastq2 (str): File with read 2 mates. Can be gzipped.
        input_bam (str): Path to an already aligend file, skipping 
            bowtie2 alignment.
        barcode (str) IUPAC barcode structure to search for in construct
            backbone. Ex: WSWSWSWSWSWSWSWSWSWS. 
            Should be specified in the guide backbone as N.
        barcode_location (int): The starting base pair of the barcode
            in the construct.
        barcode_read (int): The read in which the barcode occurs.
        guide_edit_threshold (int): The maximum number of sinlge nucleotide
            edits allowed in the guide region for an alignment.
        barcode_edit_threshold (int): The maximum number of sinlge 
            nucleotide edits allowed a barcode for an alignment.
        output_counts (str): The file to write the construct counts to.
        output_barcodes (str): The file to write the barcodes to.
        output_bam (str): The file to write the read alignments to.
        guide_length_r1 (int): The length of guide 1
        guide_length_r2 (int): The length of guide 2


    """
    def __init__(self, 
        library = None,
        fastq1 = None,
        fastq2=None,
        guide_5p_r1 = None,
        guide_3p_r1 = None,
        guide_5p_r2 = None,
        guide_3p_r2 = None,
        input_bam=None,
        barcode=None, 
        barcode_location = 175,
        barcode_read = '1',
        guide_edit_threshold=2,
        barcode_edit_threshold =0,
        guide_length_r1 = 20,
        guide_length_r2 = 20,
        sample = "Sample",
        output_counts = "/dev/stdout",
        output_bam = None,
        output_barcodes = None,
        threads = 1):


        # read in library 
        self.library = count.read_library_file(library)
        
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
        self.sample = sample 
        self.output_bam = output_bam
        self.output_counts = output_counts
        self.output_barcodes = output_barcodes

        self.threads = threads
        # get directory
        self.dir = os.getcwd()

        # are we paired?
        if self.fastq2:
            self.paired = True

        # make a temporary working directory 
        self.tmp_dir = os.path.join(self.dir,"ctg_tmp_{}".format(
            int(time.time())))
        os.makedirs(self.tmp_dir)

        self.output_summary_path = os.path.join(self.tmp_dir,
                                                "tmp.summary.csv")

        # store aligned_bam if supplied
        self.aligned_bam = input_bam

        # overwrite sample if None
        if self.sample is None:
            self.sample = os.path.splitext(os.path.basename(self.fastq1))[0]

        # check the supplied parameters
        self._check_args()

    def _check_args(self):
        """ Verify arguments pass some logic checks
        """
        if not self.fastq1:
            logging.error("Must provide at least one FASTQ file.")
            raise RuntimeError()

        if  ( not self.guide_5p_r1 ) | \
            ( not self.guide_3p_r1 ):
            logging.error("Must supply guide backbone structre: ",
                "guide_5p_r1, guide_3p_r1")
            raise RuntimeError()

        if self.fastq2:
            if  ( not self.guide_5p_r2 ) | \
                ( not self.guide_3p_r2 ):
                logging.error("Must supply guide backbone structre for read 2: ",
                    "guide_5p_r2, guide_3p_r1.")
                raise RuntimeError() # TODO

        if self.barcode:
            if (self.barcode_read == "2"):
                if (self.fastq2 is None):
                    logging.error( "Barcode was specified in read 2 ",
                                    "but no read 2 fastq was passed.")
                    raise RuntimeError()
                if  (self.barcode_location + len(self.barcode)) > \
                    (len(self.guide_5p_r2)+ \
                    len(self.guide_3p_r2)+ \
                    self.guide_length_r2):
                    logging.error( "Barcode position extends beyond ",
                                    "expected read length.")
                    raise RuntimeError()
            else:
                if (self.barcode_location +len(self.barcode)) > \
                    len(self.guide_5p_r1)+ \
                    len(self.guide_3p_r1)+ \
                    self.guide_length_r1:
                    logging.error("Barcode position extends beyond ",
                                    "expected read length.")
                    raise RuntimeError()



    def _parse_fastq_name(self):
        """ Parse fastq name to get sample, time point, and
            replicate information
            Assume fastq is formatted as:
                <SAMPLE>_[d|t]<Timepoint>_<Replicate>_...fastq
        """
        name = os.path.splitext(os.path.basename(self.fastq1))[0]
        split = name.split("_")


        if len(split)<3:
            logging.error(  "Incorrect fastq name structure, ",
                            "reverting to default. ex: SAMPLE_d3_1_r1.fastq.")
            return None

        sample = split[0]
       
        timepoint = split[1]
        if timepoint[0].upper() not in ['T','D']:
            logging.error( "Parsing fastq file name. Timepoint style string ",
                            "not in second field, reverting to default. ",
                            "Please use correct file name. ",
                            "ex: SAMPLE_d3_1_r1.fastq.")
            return None
        try:
            replicate = split[2]
            _ = int(replicate)
        except ValueError:
            logging.warning("Parsing fastq file name. Third field doesn't ",
                             "appear to be a number. Using default.")
            #return None


        fid = "_".join(split[:3])
        return fid

    def build_reference(self):
        """ Build reference from library definition.

        Writes a fastq file to and use bowtie2-build to construct kmer reference
        to align sequencing reads against.

        Uses object level parameters
        library (dict): dictionary representation of the dual CRISPR guides 
        see: count.read_library_file
        g_5p_r1 (str): the 5' sequence of the construct upstream of the 
        first guide.
        g_3p_r1 (str): the 3' sequence of the construct downstream of the
        first guide.
        g_5p_r2 (str): the 5' sequence of the construct upstream of the 
        second guide.
        g_3p_r2 (str): the 3' sequence of the construct downstream of the 
        second guide.
        reference_fasta (str): the name of the reference file
        tmp_dir (str): directory path to write the reference file to

        Creates a bowtie2 index file.
        """
        if self.aligned_bam:
            logging.warning("An aligned bam was provided as an input, ",
                            "negating the need to build a reference.")
            return 0

        self.reference = os.path.join(self.tmp_dir, "guide_reference.fa")

        # parse library definition and write out guide refernece
        count.build_guide_reference(  self.library,
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
        """ 
        Align single reads containing a guide sequence with bowtie2 to the 
        guide based reference.  Additionally extracting the guide sequence, 
        calculating the edit distance in this guide and add both as tags to
        bam file containing the alignments.

        Optianlly, extract a barcode sequence and calculate the edit 
            distance, again adding both as bam tags.

        Optionally, add bitwise flags to the alingments.

        Args:
            reference (str): Path to a reference fasta file to build 
                bowtie2 index from
            bt2_index_base (str): Base path containing the bowtie2 index
            guide_start (int): Position from the 5' end of the read the 
                guide sequence starts in the reads 
            guide_length (int): The expected length of the guide sequence
            expected_barcode: The expected barcode sequence in IUPAC format. 
                ex: WSWS for [A|T][G|C][A|T][G|C][A|T][G|C]
                so AGTC would be extracted
            barcode_start (int): Position from the 5' end of the read the 
                barcode sequence starts in the reads 
            binary_flags_to_add (int): Bitwise flags to add to the reads. 
                ex: 65 for a paired read 1
            bam (str): the output bam file to write the single end reads to
    
        Returns:
            aligned_bam (path to file): Creates a bam file containing 
                aligned reads with probe assignments.
        """
        if self.aligned_bam:
            logging.warning("An aligned bam was provided as an input, ",
                            "negating the need to align.")
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
        ret_code = align.bowtie2_align(self.fastq1, 
            self.reference, 
            "--quiet", 
            "--very-sensitive",
            bam=bam1,
            binary_flags_to_add=r1_flags_to_add,
            guide_start=len(self.guide_5p_r1), 
            guide_length=self.guide_length_r1,
            expected_barcode=barcode,
            barcode_start=barcode_location,
            threads=self.threads )
        if ret_code != 0:
            logging.error("Bowtie2 error. See trace.")
            self.cleanup()
            sys.exit(ret_code)

        if self.fastq2:
            if self.barcode_read == "2":
                barcode=self.barcode
                barcode_location = self.barcode_location
            else:
                barcode = None
                barcode_location = None
            # align read 2 
            bam2 = os.path.join(self.tmp_dir, "tmp.r2.bam")
            align.bowtie2_align(self.fastq2, 
                self.reference, 
                "--quiet", 
                "--very-sensitive",
                bam=bam2,
                binary_flags_to_add=r2_flags_to_add,
                guide_start=len(self.guide_5p_r2), 
                guide_length=self.guide_length_r2,
                expected_barcode=barcode,
                barcode_start=barcode_location,
                threads=self.threads)

            # merge
            self.aligned_bam = os.path.join(self.tmp_dir, "tmp.bam")
            ret_code = align.merge_fixmate(bam1,
                                bam2,
                                self.aligned_bam,
                                tmp=os.path.join(self.tmp_dir,"sort"))
            if ret_code != 0:
                logging.error("Bowtie2 error. See trace.")
                self.cleanup()
                sys.exit(ret_code)
            # remove temporary read level bams
            os.remove(bam1)
            os.remove(bam2)

        else:
            # if no read 2, just keep read 1 aligned bam and store path
            self.aligned_bam = bam1

    def count(self):
        """ Count construct pairs in aligned bam file.

        Uses object level parameters:
            bam_file_path (str): the bam file path containing the reads to 
                count
            library (dict): dictionary representation of the dual CRISPR 
                guides see: count.read_library_file
            guide_edit_threshold (int): the maximum number of single base pair
                edits to allow in the guide region of aligned reads.
            output_counts_path (str): the path to write the construct counts 
                to. Default is standard out.

        Creates a construct counts file.
        """
        log = logging.getLogger()
        if self.aligned_bam is None:
            raise RuntimeError("Must align reads before counting.")

        log.debug("Aligned bam: {}".format(self.aligned_bam))
        
        if (not os.path.exists(self.aligned_bam)):
            raise RuntimeError("Bam does not exists. Must first align fastqs.")

        log.info("Couting aligned constructs.")
        if self.barcode:
            count.count_good_constructs_and_barcodes(self.aligned_bam,
                library = self.library, 
                guide_edit_threshold = self.guide_edit_threshold,
                barcode_edit_threshold = self.barcode_edit_threshold,
                sample = self.sample,
                output_counts_path = self.output_counts,
                output_barcodes_path = self.output_barcodes)
        else:
            count.count_good_constructs(self.aligned_bam,
                library = self.library, 
                sample = self.sample,
                guide_edit_threshold = self.guide_edit_threshold,
                output_counts_path = self.output_counts)


        return

    def cleanup(self):
        """ Remove tmp directory created during construct counting.

        Uses object level parameters:
            output_bam (str): The file to store the read alignments in.
            aligned_bam (str): The file alignments were written to.
            tmp_dir (str): The temporary working directory used.
        """
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
