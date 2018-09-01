import sys
import re
import pysam
import shutil
import os 
import gzip 

from multiprocessing import Pool
from optparse import OptionParser
from collections import Counter
from contextlib import contextmanager

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to process bam files split by chromosome and add a new tag based on the haplotype"
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the .bam file to parse")
opts.add_option("--wl1", help="Filepath to first haplotype whitelist")
opts.add_option("--wl2", help="Filepath to second haplotype whitelist")
opts.add_option("--bam-tag", "-b", help="Two characters to be associated with the bam tag that will be appended")

opts.add_option("--out", "-o", help="Annotated bam file output path")

options, arguments = opts.parse_args()

inputbamname = options.input
outputbamname = options.out
wl1_file = options.wl1
wl2_file = options.wl2
bam_tag = options.bam_tag

# Read in the whitelists
with open(wl1_file) as wl1_file_h:
	content = wl1_file_h.readlines()
	wl1 = [x.strip() for x in content] 

with open(wl2_file) as wl2_file_h:
	content = wl2_file_h.readlines()
	wl2 = [x.strip() for x in content] 

# BAM I/O
bam = pysam.AlignmentFile(inputbamname, "rb")
out = pysam.AlignmentFile(outputbamname, "wb", template = bam)

# Loop over bam and extract the sequence 
for read in bam:
	read_name = read.query_name
	# If read barcode is in whitelist, then write it out
	tag = 0
	if(read_name in wl1):
		tag = 1
	elif(read_name in wl2):
		tag = 2
	else:
		tag = 0 
	read.tags = read.tags + [(bam_tag, tag)]
	out.write(read)
	
bam.close()
out.close()