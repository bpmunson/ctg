"""
Wrapper for bowtie2 functions
"""

import os
import subprocess
import logging


import regex as re # use for fuzzy matching 
import numpy as np






# Command line functions
def bowtie2_build(reference_in, bt2_index_base, *args, **kwargs):
	"""
	reference_in - reference fastq to index
	bt2_index_base - optional prefix to store index
	args - addtional bowtie2-build flags
	kwargs - additional bowtie2-build keword args. examples: o=5, ftabchars=10

	No validation of agruments is performed.
	"""

	try:
		subprocess.check_output('bowtie2-build --version', shell=True)
	except OSError:
		raise RuntimeError('Commnad bowtie2-build not found. Please install bowtie2.')

	call = ["bowtie2-build"]
	# add optional arguments
	for arg in args:
		call.append(arg)
	for k, v in kwargs.items():
		if len(k)==1: # single 
			call.append("-{} {}".format(k, v))
		else:
			call.append("--{} {}".format(k, v))
	if bt2_index_base:
		call.append(bt2_index_base)
	call.append(reference_in)
	#print("call: ", call)
	logging.info("Building Bowtie2 Index.")
	subprocess.check_call(call)


	return 0

def bowtie2_align(fastq, bt2_index_base, *args, 
	guide_start=None, guide_length=None,
	expected_barcode=None, barcode_start=None,
	binary_flags_to_add = None, bam=None, **kwargs):
	""" Align single reads with bowtie2 to reference

		if qsub is specified we instead just write out shells for the user to submit manually
	"""
	try:
		subprocess.check_output('bowtie2 --version', shell=True)
	except OSError:
		raise RuntimeError('Commnad bowtie2 not found. Please install bowtie2.')

	call = ['bowtie2']
	# add optional arguments
	for arg in args:
		call.append("{}".format(arg))
	for k, v in kwargs.items():
		call.append("--{} {}".format(k, v))

	call.append("-x {}".format(bt2_index_base))
	call.append("-U {}".format(fastq))

	# add 
	if guide_start or expected_barcode:
		call.append("| add_tags.py")
		if guide_start:
			if guide_length is None:
				raise RuntimeError()
			call.append("--guide_start {} --guide_length {}".format(guide_start, guide_length))
		if expected_barcode:
			if barcode_start is None:
				raise RuntimeError()
			call.append("--expected_barcode {} --barcode_start {}".format(expected_barcode, barcode_start))
		if binary_flags_to_add:
			call.append("--flag {}".format(binary_flags_to_add))


	# add conversion to bam format
	if bam is None:
		bam = "{}.bam".format(os.path.splitext(fastq)[0])

	call.append("| samtools view -Sb - > {}".format(bam))

	# merge call and execute	
	call = " ".join(call)
	logging.info("Running Bowtie2 Alignment: {}".format(call))
	subprocess.check_output(call, shell=True)

	return 0

def merge_fixmate(r1_bam, r2_bam, out_bam, tmp="./sort"):
	""" Merge the read 1 and 2 bams, queryname sort, and run fixmate
	"""
	try:
		subprocess.check_output('samtools --version', shell=True)
	except OSError:
		raise RuntimeError('Commnad bowtie2 not found. Please install samtools.')	
	call = ["samtools cat", r1_bam, r2_bam, "| samtools sort", "-T", tmp, "-n - | samtools fixmate -", out_bam]

	# merge call and execute	
	call = " ".join(call)
	logging.info("Merging and Fixing Mate Information: {}".format(call))
	subprocess.check_output(call, shell=True)
	return 0







# def bowtie2_align_old(fastq, bt2_index_base, *args, bam=None, binary_flags_to_add = None, **kwargs):
# 	""" Align single reads with bowtie2 to reference

# 		if qsub is specified we instead just write out shells for the user to submit manually
# 	"""
# 	try:
# 		subprocess.check_output('bowtie2 --version', shell=True)
# 	except OSError:
# 		raise RuntimeError('Commnad bowtie2 not found. Please install bowtie2.')

# 	call = ['bowtie2']
# 	# add optional arguments
# 	for arg in args:
# 		call.append("{}".format(arg))
# 	for k, v in kwargs.items():
# 		call.append("--{} {}".format(k, v))

# 	call.append("-x {}".format(bt2_index_base))
# 	call.append("-U {}".format(fastq))
# 	# add awk call to change binary flags 
# 	if binary_flags_to_add:
# 		call.append("| awk '{{if($1 !~ /^@/){{$2+={0}}};print $0}}' OFS='\\t'".format(binary_flags_to_add))
# 	# add conversion to bam format
# 	if bam is None:
# 		bam = "{}.bam".format(os.path.splitext(fastq)[0])

# 	call.append("| samtools view -Sb - > {}".format(bam))

# 	# merge call and execute	
# 	call = " ".join(call)
# 	logging.info("Running Bowtie2 Alignment: {}".format(call))
# 	subprocess.check_output(call, shell=True)
#
#	return 0