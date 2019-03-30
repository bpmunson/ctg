__all__ = ['core','Scorer','Counter']

#from ctg import *
import configargparse
import sys
import os
import logging
import pysam

from ctg.Counter import Counter
from ctg.Scorer import Scorer
from ctg.core.aggregate import aggregate_counts
from ctg.core.count import add_tags
from ctg.version import __version__

# Parsers
def count_parser(subparser):
    count_subparser = subparser.add_parser('count', help="Count construct reads from a fastq file.")


    # guide structure
    count_subparser_req_read=count_subparser.add_argument_group(title='Required arguments for read structure',description='')
    count_subparser_req_read.add_argument("--guide_5p_r1", action="store", required=True, type=str, help="Expected 5' end of read 1.")
    count_subparser_req_read.add_argument("--guide_3p_r1", action="store", required=True, type=str, help="Expected 3' end of read 1.")

    count_subparser_req=count_subparser.add_argument_group(title='Required arguments',description='')
    count_subparser_req.add_argument("-l","--library", action="store", required=True, help='Library defintion file in csv format.  ')
    count_subparser_req_me = count_subparser_req.add_mutually_exclusive_group(required=True)
    count_subparser_req_me.add_argument("-1","--fastq1", action="store", default=None, help='File with read 1 mates.  Can be gzipped.')





    count_subparser_opt_rt=count_subparser.add_argument_group(title='Optional run time arguments',description='')
    count_subparser_opt_rt.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose output.")
    count_subparser_opt_rt.add_argument("-q", "--quiet", action="store_true", default=False, help="Supress all warnings and info.  Supersedes --verbose.")
    count_subparser_opt_rt.add_argument("-d", "--debug", action="store_true", default=False, help="Debug output. Supersedes --quiet.")
    count_subparser_opt_rt.add_argument("--threads", action="store", required=False, default=1, help='Number of threads to use.')


    # required arguments
    # fastqs 

    # optional arguments
    count_subparser_opt=count_subparser.add_argument_group(title='Optional input arguments',description='')
    count_subparser_opt.add_argument("--sample", action="store", required=False, default=None, help="Sample name to add to counts file.")
    count_subparser_opt.add_argument("-c", "--config", required=False, is_config_file=True, help="Path to config file.")


    # barcode structure : # TODO make these all required if one is set
    count_subparser_opt_barcode=count_subparser.add_argument_group(title='Optional arguments for finding barcodes/UMI',description='')

    count_subparser_opt_barcode.add_argument("--barcode", action="store", default=None, help="IUPAC barcode structure to search for in construct backbone. Ex: WSWSWSWSWSWSWSWSWSWS. Should be specified in the guide backbone as N.")
    count_subparser_opt_barcode.add_argument("--barcode_location", action="store", default=175, type=int, help="The starting base pair of the barcode in the construct.")
    count_subparser_opt_barcode.add_argument("--barcode_read", action="store", default='1', type=str, help="The read in which the barcode occurs.")


    # Threshold for counting
    count_subparser_opt_count=count_subparser.add_argument_group(title='Optional arguments for counting constructs',description='')
    count_subparser_opt_count.add_argument("--guide_edit_threshold", action="store", default=2, type=int, help="The maximum number of sinlge nucleotide edits allowed in the guide region for an alignment." )
    count_subparser_opt_count.add_argument("--barcode_edit_threshold", action="store", default=2, type=int, help="The maximum number of sinlge nucleotide edits allowed a barcode for an alignment." )

    # output files
    count_subparser_opt_output=count_subparser.add_argument_group(title='Optional arguments for counting outputs',description='')
    count_subparser_opt_output.add_argument("--output_counts", action="store", default="/dev/stdout", help="The file to write the construct counts to." )
    count_subparser_opt_output.add_argument("--output_barcodes", action="store", default=None, help="The file to write the barcodes to." )
    count_subparser_opt_output.add_argument("--output_bam", action="store", default=None, help="The file to write the read alignments to." )
    count_subparser_opt_output.add_argument("--leave_tmp", action="store_true", default=False, help="Do not clean up the tmp working directory.")

    count_subparser_opt.add_argument("--guide_5p_r2", action="store", required=False, type=str, help="Expected 5' end of read 2.")
    count_subparser_opt.add_argument("--guide_3p_r2", action="store", required=False, type=str, help="Expected 3' end of read 2.")
    count_subparser_opt.add_argument("-2","--fastq2", action="store", required=False, help='File with read 2 mates.  Can be gzipped.')
    count_subparser_opt.add_argument("--input_bam", action="store", default=None, help="Path to an already aligend file, skipping bowtie2 alignment.")

def add_tags_parser(subparser):
    """ subparser for add tags
    """
    add_tags_subparser = subparser.add_parser("add_tags", help="Add tags to a bam file")
    add_tags_subparser.add_argument("--flag", action="store", default=None, type=int, help="Bitwise flags to add to all reads.")
    add_tags_subparser.add_argument("--expected_barcode", action="store", type=str, default=None,required=False, help="Position of the guide in the read.")
    add_tags_subparser.add_argument("--barcode_start",action="store",type=int,default=None,required=False,help="Position of barcode in read")
    add_tags_subparser.add_argument("--guide_start", action="store", type=int,help="Position of the guide in the read.")
    add_tags_subparser.add_argument("--guide_length", action="store",type=int, help="Position of the guide in the read.")

def aggregate_parser(subparser):
    """ Add options for aggregate function
    """
    aggregate_subparser = subparser.add_parser('aggregate', help="Aggregate counts files.")
    aggregate_subparser_rt=aggregate_subparser.add_argument_group(title='Optional run time arguments',description='')
    aggregate_subparser_rt.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose output.")
    aggregate_subparser_rt.add_argument("-q", "--quiet", action="store_true", default=False, help="Supress all warnings and info.  Supersedes --verbose.")
    aggregate_subparser_rt.add_argument("-d", "--debug", action="store_true", default=False, help="Debug output. Supersedes --quiet.")
    aggregate_subparser_rt.add_argument("--threads", action="store", required=False, default=1, help='Number of threads to use.')

    # required arguments
    # fastqs 
    aggregate_subparser_req=aggregate_subparser.add_argument_group(title='Required arguments',description='')

    # library defintion
    #aggregate_subparser_req.add_argument("--library", action="store", required=True, help='Library defintion file in csv format.  ')
    aggregate_subparser_req.add_argument("--counts_files", action="store", required=True, nargs="*", help='Space separated list of counts files created by ctg count.')

    # optional Arguments
    aggregate_subparser_opt=aggregate_subparser.add_argument_group(title='Optional input arguments',description='')
    aggregate_subparser_opt.add_argument("--names", action="store", default=None, nargs="+", help="Space separated list of names to use as column headers in output.")
    aggregate_subparser_opt.add_argument("--output", action="store", default="/dev/stdout", help="Output path")

def score_parser(subparser):
    score_parser = subparser.add_parser('score', help="Calculate genetic interaction scores from construct counts.")

    score_parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose output.")
    score_parser.add_argument("-q", "--quiet", action="store_true", default=False, help="Supress all warnings and info.  Supersedes --verbose.")
    score_parser.add_argument("-d", "--debug", action="store_true", default=False, help="Debug output. Supersedes --quiet.")
    score_parser.add_argument("--threads", action="store", required=False, default=1, help='Number of threads to use.')
    score_parser.add_argument("--testing", action="store_true", default=False, help="Flag for testing purposes only.")
    # fit_ac_fc
    score_parser.add_argument("--min_time_points", action="store", default=2, type=int, help="Minimum number of timepoints to use in fitness estimation.")
    # irls
    score_parser.add_argument("--bi_weight_stds", action="store", default=2, type=float, help="Number of standard deviation to use in Tukey biweight normalization.")
    score_parser.add_argument("--tol", action="store", default=1e-3, type=float, help="Relative error tolerance")
    score_parser.add_argument("--max_irls_iter", action="store", default=50, type=int, help="Maximum IRLS iterations to perform")
    # weighted pi
    score_parser.add_argument("--n_probes_per_target", action="store", default=2, type=int, help="Maximum number of probes per target to use.") 
    score_parser.add_argument("--iterations", action="store", type=int, default=2, help="Number of bootstrap iterations to perform")
    score_parser.add_argument("--null_target_id", action="store", default=None, help="Target/Gene name which corresponds to the null target")
    # input
    score_parser.add_argument("-c", "--time_point_counts", action="store", default=None, required=True, help="Path to timepoint counts file.")
    score_parser.add_argument("-a", "--abundance_file", action="store", default=None, required=False, help="Comma separated list of timepoints to use.")
    # output
    score_parser.add_argument("-o", "--output", action="store", default=None, help="Output results path")
    score_parser.add_argument("-p", "--pickle_output", action="store", default=None, help="Outuput a pickled object to this path.")
    
    # extra 
    #score_parser.add_argument("--dont_use_full_dataset_for_ranking", action="store_true", default=False,help="Flag to use all the dataset for probe rankings in bootstraps.")
    #score_parser.add_argument("-a", "--abundances", action="store", default=None, required=False, help="Path to timepoint counts file.")
    #score_parser.add_argument("--only_use_sampled_fc_for_pi_estimates", action="store_true", default=False, help="Don't estimate pi scores for all the constructs in bootstraps")
    #score_parser.add_argument("--make_single_gene_screen", action="store_true", default=False,  help="Filter out everything but gene by nontargeting")
    #score_parser.add_argument("--null_aware", action="store_true", default=True, help="Treat the null probes different in rankings.")

def ctg_parseargs():
    parser = configargparse.ArgumentParser(description="CTG analysis.")

    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--verbose", action="store_true", default=False, help="Verbose output.")
    parser.add_argument("--quiet", action="store_true", default=False, help="Supress all warnings and info.  Superceeds --verbose.")
    parser.add_argument("--debug", action="store_true", default=False, help="Debug output. ")


    subparser = parser.add_subparsers(help="Commands Available", dest="command")

    # add subparser arguments
    count_parser(subparser)
    aggregate_parser(subparser)
    score_parser(subparser)
    add_tags_parser(subparser)


    #args = parser.parse_args()

    return parser

# Mains
def counting_main(options):
    """ Main to run the counting pipeline
    """
    # build bowtie refernce from library defintion
    counter = Counter(
            library = options.library, 
            fastq1 = options.fastq1,
            input_bam = options.input_bam,
            fastq2 = options.fastq2,
            barcode = options.barcode,
            barcode_location = options.barcode_location,
            barcode_read = options.barcode_read,
            guide_edit_threshold = options.guide_edit_threshold,
            barcode_edit_threshold = options.barcode_edit_threshold,
            guide_5p_r1 = options.guide_5p_r1,
            guide_3p_r1 = options.guide_3p_r1,
            guide_5p_r2 = options.guide_5p_r2,
            guide_3p_r2 = options.guide_3p_r2,
            sample = options.sample,
            output_counts = options.output_counts,
            output_bam = options.output_bam,
            output_barcodes = options.output_barcodes,
            leave_tmp = options.leave_tmp,
            threads= options.threads)

    # run pipeline
    counter.build_reference()
    counter.align_reads()
    counter.count()
    if not options.leave_tmp:
        counter.cleanup()
    return 0

def aggregate_main(options):
    """ Main to aggregate read counts 
    """
    #print(options)
    #count.aggregate_counts(options.library, options.counts_files, options.output, options.names )
    aggregate_counts(options.counts_files, output_file = options.output, sample_names = options.names)
    return 0

def scoring_main(options):
    """ Main to run the scoring pipeline
    """

    scorer = Scorer(options.time_point_counts, 
        verbose = options.verbose,
        min_time_points = options.min_time_points,
        bi_weight_stds = options.bi_weight_stds,
        tol = options.tol,
        maxiter = options.max_irls_iter,
        n_probes_per_target = options.n_probes_per_target,
        null_target_id = options.null_target_id,
        niter = options.iterations,
        testing = options.testing,
        output = options.output,
        pickle_output = options.pickle_output,
        threads = options.threads)
        #use_full_dataset_for_ranking = not options.dont_use_full_dataset_for_ranking,

    # run pipeline
    scorer.run_construct_fitting()
    scorer.run_pi_score_calculation()
    scorer.run_weighting()
    scorer.run_sampling()
    scorer.summarize()
    if scorer.pickle_output is not None:
        scorer.pickle()

    return 0

def add_tags_main(options):
    """ Main for add_tags function
    """

    bam = pysam.AlignmentFile("-", "r")

    # get header
    header= bam.header.to_dict()
    new_pg = {'CL':' '.join(sys.argv), 'ID':'CTG', 'PN':'add_tags', 'VN':__version__}
    header['PG'].append(new_pg)
    output = pysam.AlignmentFile("-", "w", header=header)

    while True:
        logging.info("Reading in bam file and adding tags.")
        try:
            read = next(bam)
        except StopIteration:
            break

        if not read.is_unmapped:
            try:
                read = add_tags(read, options.guide_start, options.guide_length, options.expected_barcode, options.barcode_start, options.flag)
            except IndexError:
                logging.error("Error in Read: {}".format(read.qname))
                read = add_tags(read, flag=options.flag)
                output.write(read)
                exit()
        else:
            read = add_tags(read, flag=options.flag)

        try:
            output.write(read)
        except OSError:
            break

    logging.info("Finished")
    bam.close()
    output.close()

def main():
    """ Main
    """
    parser = ctg_parseargs()
    args = parser.parse_args()
    
    # set up logger
    level = logging.WARNING
    if args.quiet:
        level = logging.ERROR
    if args.verbose:
        level = logging.INFO
    if args.debug:
        level = logging.DEBUG 


    logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=level) 
    logger = logging.getLogger()

    if args.command == None:
        # no command specified
        parser.print_help()
        sys.exit(0)

    if args.command == "count":
        logger.info('Startng CTG counting pipeline.')
        counting_main(args)
        logger.info("Exiting")
        sys.exit(0)

    if args.command == "score":
        logger.info('Starting CTG scoring pipeline.')
        scoring_main(args)
        logger.info("Exiting")
        sys.exit(0)

    if args.command == "aggregate":
        logger.info("Starting CTG aggregation.")
        aggregate_main(args)
        logger.info("Exiting")
        sys.exit(0)

    if args.command == "add_tags":
        logger.info("Adding bam tags.")
        add_tags_main(args)
        sys.exit(0)

    return 0 

if __name__ == '__main__':
    main()


