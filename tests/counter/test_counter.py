import os
import logging
import numpy as np
import pandas as pd
import ctg
import filecmp

def _compare_to_benchmark(testing, reference):
    """ Compare a numpy array to the benchmark """
    assert np.allclose(testing, reference)

def _compare_files(testing, reference):
	""" Ensure two files are identical in content """
	res = filecmp.cmp(testing, reference)
	assert(res)


def _compare_counts(counter):
	""" Compare core counts """	
	pass

def main():
	""" Run counting and perform tests """
	files = {
		"counts": "benchmarks/testing_counts_benchmark.csv",
	}

	# set up counter args
	config = ctg.core.config.config
	fastq1 = os.path.join(config.A549_test, "A549-CV4-1000_d21_1_R1.fastq")
	fastq2 = os.path.join(config.A549_test, "A549-CV4-1000_d21_1_R2.fastq")
	library = os.path.join(config.A549_test, "library_definition.CV4.csv")
	guide_5p_r1 = "TATATATCTTGTGGAAAGGACGAAACACCG"
	guide_3p_r1 = "GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGGCCACCTTAACACGCGATGATATTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCTATTACGAGCGCTTGGAT"
	guide_5p_r2 = "CTTGGAGAAAAGCCTTGTTTG"
	guide_3p_r2 = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG"
	guide_edit_threshold=2
	output_counts = "A549-CV4-d21_1_counts.txt"


	# initialize counter
	counter = ctg.Counter(library, fastq1, fastq2=fastq2,
		guide_5p_r1=guide_5p_r1,
		guide_5p_r2=guide_5p_r2,
		guide_3p_r1=guide_3p_r1,
		guide_3p_r2=guide_3p_r2,
		guide_edit_threshold=guide_edit_threshold,
		output_counts = output_counts)

	# run counting pipeline
	counter.build_reference()
	counter.align_reads()
	counter.count()

	# run comparisons
	_compare_files(output_counts, files['counts'])

	# remove  files
	counter.cleanup()
	os.remove(output_counts)
	return 0


if __name__ == "__main__":
    # set up logger
    logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s',
                        level=logging.INFO) 
    # Run Tests
    ret_code = main()

    if ret_code == 0:
        logging.info("All Tests Passed")
    else:
    	logging.error("Test(s) Failed")
