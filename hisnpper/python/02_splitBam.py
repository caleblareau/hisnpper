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
usage = "usage: %prog [options] [inputs] Basic script to split bam files by chromosome name and sort by name."
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the .bam file to parse")
opts.add_option("--out", "-o", help="Path to the output directory for these files")
opts.add_option("--ncores", default = 4, help="Number of cores for parallel processing.")
opts.add_option("--chrfile", help="Path to chrs")
opts.add_option("--barcode-tag", default = 'XB', help="Name of the first .bam file")


options, arguments = opts.parse_args()

bamname = options.input
out = options.out
barcode_tag = str(options.barcode_tag)
cpu = int(options.ncores)
chrs = [line.strip() for line in open(options.chrfile, 'r')]

def getBarcode(intags):
	'''
	Parse out the barcode per-read
	'''
	if(barcode_tag == "none"):
		return("none")
	
	for tg in intags:
		if(barcode_tag == tg[0]):
			return(tg[1])
	return("NA")


# Handle temporary directory structure
temp = out + "/temp"
temp_split = temp + "/01_split"

#---------------------------------------------------------
# Function for writing the read name and the bead barcode ID for later use
#---------------------------------------------------------
def writeBamRead(two):
	chr = two[0]
	filename = two[1]
	
	# Iterate through bam file
	bam = pysam.AlignmentFile(bamname,'rb')
	Itr = bam.fetch(str(chr),multiple_iterators=True)
	with gzip.open(filename, 'wt') as out_write:
		for read in Itr:
		
			# Write reads + barcodes
			read_barcode = getBarcode(read.tags)
			read_name = read.query_name
			value = str(read_name) + "\t" + str(read_barcode) + "\n"
			out_write.write(value)
	bam.close()
	
	# Split into per-chromosome bam files
	new_bam = filename.replace(".read_barcode.tsv.gz", ".bam")
	
	# Execute a split via samtools
	split_cmd = "samtools view -b " + bamname + " " + chr + " > " + new_bam
	os.system(split_cmd)
	pysam.index(new_bam)

read_barcode_file = [temp_split + "/" + "splitBam" + "." + chr + ".read_barcode." +"tsv.gz" for chr in chrs]

pool = Pool(processes=cpu)
toy_out = pool.map(writeBamRead, zip(chrs, read_barcode_file))
pool.close()




