
import os
import time
import align
import shutil
import logging
import count

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

        log.info("Couting aligned constructs.")
        count.count_good_constructs(self.aligned_bam,
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
