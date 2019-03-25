import sys
import re
import pysam
import shutil
import os 
import gzip 

from optparse import OptionParser
from collections import Counter
from contextlib import contextmanager
from string import digits

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to process raw data from bam split."
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the data file to parse")
opts.add_option("--output", "-o", help="Annotated bam file output path")
opts.add_option("--cutoff", "-c", help="Cutoff for alignment quality")
opts.add_option("--complement", "-r", help="Compliment the bases (both ref and alt)-- applicable for - strand")

options, arguments = opts.parse_args()

input_file_name = options.input
output_file_name = options.output
cutoff = float(options.cutoff)
complement = options.complement == "yes"

# Reformat the calmd raw stuff

remove_digits = str.maketrans('', '', digits)
complement_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N' : 'N'}

# Function to offset and mask read
def process_cigar_soft_mask(read, cigar, pos):
	
	# Case 1: less than 10 bp soft clipping to start
	if(cigar[1] == "S"):
		nbp = int(cigar[0])
		read = "="*nbp + read[nbp:]
		pos = pos - nbp
		
	# Case 2: at least 10 bp soft clipping to start read
	if(cigar[2] == "S"):
		nbp = int(cigar[0:2])
		read = "="*nbp + read[nbp:]
		pos = pos - nbp
		
	# Case 3: end of read (cases are not mutually exclusive)
	if(cigar[-1] == "S"):

		nbp = 0
		# Facet a-- at least 10 bp
		if(cigar[-3] in digits):
			nbp = int(cigar[-3:-1])
		else:
			nbp = int(cigar[-2])
		read = read[:(-1*nbp)] + "="*nbp
		# No need to offset the read since reading from the left
	return(read, pos)
	

with open(input_file_name) as ino:
	with open(output_file_name, 'w') as ono:
		for line in ino:
			v = line.split("\t")
			if(len(v) > 10): # filters out header lines
				align_quality = float(v[4])
				read = v[9]
				cigar = v[5]
				pos = int(v[3])
				
				# Handle soft-masking by overwriting soft-masked bases with "=" and shifting index, if needed
				if("S" in cigar):
					read, pos = process_cigar_soft_mask(read, cigar, pos)

				# Only process read if it passes alignment quality and has a mismatch with the reference
				if((align_quality > cutoff) & (read.count("=") < len(read)) & ("I" not in cigar) & ("D" not in cigar)):
					
					# Most of the try / catching occurs to facilitate handling reads with indels
					try:
						read_name = v[0]
						chromosome = v[2]
						
						BQ = v[10]

						# Pull the MD tag -- note the -7 is hard coded to look in the final seven elements
						# This is to pull the reference genome location
						MD_tag = (list(filter(lambda x:'MD:' in x, v[-7:]))[0])[5:]
						ref_ordered = MD_tag.translate(remove_digits).strip()
			
						# Make two indices; i indexes which mismatch we are in; j indexes where in the read we are
						i = 0; j = 0
						for c in read:
							if c != "=":
						
								# Make new attribtues for single base pair
								ref = ref_ordered[i]
								alt = read[j]
								if(complement):
									ref = complement_dict[ref]
									alt = complement_dict[alt]
							
								# Prepare final output	
								edit = ref + "_" + alt
								pos_new = str(int(pos) +j)
								BQ_out = BQ[j]
						
								# Export the result
								ono.write(read_name + "\t" + chromosome + "\t" + pos_new + "\t" + edit + "\t" + BQ_out + "\n")
		
								i += 1 # increase 1 which mutation we are looking at
							j += 1 # walk forward 1 position in read 
					except:
						pass

