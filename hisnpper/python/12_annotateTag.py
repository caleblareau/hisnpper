import sys
import re
import pysam
import shutil
import os 
import gzip 

from optparse import OptionParser
from collections import Counter
from contextlib import contextmanager

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to process bam files split by chromosome and add a new tag based on the haplotype"
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the .bam file to parse")
opts.add_option("--assign-table", "-t", help="Filepath to the table with haplotype assignments")
opts.add_option("--bam-tag", "-b", help="Two characters to be associated with the bam tag that will be appended")
opts.add_option("--out", "-o", help="Annotated bam file output path")
opts.add_option("--cutoff", "-c", help="Cutoff for haplotype assignment")

options, arguments = opts.parse_args()

inputbamname = options.input
outputbamname = options.out
at_file = options.assign_table
bam_tag = options.bam_tag
cutoff = float(options.cutoff)

# Read in the whitelists
d = {}
with open(at_file) as file_h:
	for line in file_h:
		(read, haplotype, value) = line.split()
		if(float(value) > cutoff):
			d[read] = haplotype
 
def get_haplotype(read_name):
	'''
	Parse out the bead barcode per-read
	'''
	if read_name in d:
		return(d[read_name])
	else:
		return("NA")

# BAM I/O
bam = pysam.AlignmentFile(inputbamname, "rb")
out = pysam.AlignmentFile(outputbamname, "wb", template = bam)

# Loop over bam and extract the sequence 
for read in bam:
	read_name = read.query_name
	
	# If read barcode has been assigned, then write it out
	tag = get_haplotype(read_name)
	read.tags = read.tags + [(bam_tag, tag)]
	out.write(read)
	
bam.close()
out.close()