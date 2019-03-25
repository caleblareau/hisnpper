import click
import os
import os.path
import sys
import shutil
import yaml
import random
import string
import itertools
import time
import pysam
import glob
import csv
import re
from itertools import groupby

from pkg_resources import get_distribution
from subprocess import call, check_call
from .frombapHelp import *
from .hisnpperProjectClass import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs

@click.command()
@click.version_option()

@click.option('--bamfile', '-b',required = True, help='Bam file that will be parsed / annotated.')
@click.option('--fasta', '-f', required = True, help='Input fasta file; needs to be hard-masked at relevant base pairs')
@click.option('--name', '-n', default="default", help='Name for all of the output files (default: uses the .bam prefix)')
@click.option('--output', '-o', default="hisnpper_out", help='Output directory for analysis; see documentation.')

@click.option('--edit', '-e', required = True, help='Type of nucleotide modification to search for. Should be in the format of REF_ALT (e.g. C_T).')
@click.option('--min-aq', '-ma', default=20, help='Minimum alignment quality for read to be considered.')

@click.option('--barcode-tag', '-bt', default="none", help='Two letter tag that indicates the single-cell ID; by default, "none". This is essential for single-cell analyses.')

@click.option('--keep-positions', '-k', default = "none", help='Two-column (chr\tbp) file of positions to be considered, exclusively.')
@click.option('--remove-positions', '-r', default = "none", help='Two-column (chr\tbp) file of positions to be removed.')

@click.option('--ncores', '-c', default=2, help='Number of cores to be used in analysis')
@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')
@click.option('--snake-stdout', '-so', is_flag=True, help='Write snakemake log to sdout rather than a file')

def main(bamfile, fasta, name, output,
	edit, min_aq, barcode_tag, keep_positions, remove_positions,
	ncores, keep_temp_files, snake_stdout):
	
	"""
	hisnpper-edit: deriving nucleotide-specific inferences from sequencing data.\n
	Caleb Lareau, clareau <at> broadinstitute <dot> org \n
	"""
	
	__version__ = get_distribution('hisnpper').version
	script_dir = os.path.dirname(os.path.realpath(__file__))
	click.echo(gettime() + "Starting hisnpper-edit v%s" % __version__)


	# Determine chromosomes in bam file
	bam_chrs = []
	for l in pysam.idxstats(bamfile).split('\n'):
		t = l.split("\t")
		if(len(t) > 3):
			if(float(t[2]) > 0):
				bam_chrs.append(t[0])

	# Make the project
	p = hisnpperEditProject(script_dir,
		bamfile, fasta, name, output,
		edit, min_aq, barcode_tag, keep_positions, remove_positions,
		ncores, keep_temp_files)

	# Make output folders
	of = output; fin = of; temp = of + "/temp"; logs = of + "/logs";
	temp_split = temp + "/01_split"; temp_namesort = temp + "/02_frombam"
	temp_processed = temp + "/03_processed"
	
	folders = [of, fin, temp, logs, of + "/.internal", logs + "/samtools", logs + "/readstats",
			temp_split, temp_namesort, temp_processed, 
			of + "/.internal/parseltongue", of + "/.internal/samples"]
	mkfolderout = [make_folder(x) for x in folders]

	# Make internal README files
	if not os.path.exists(of + "/.internal/README"):
		with open(of + "/.internal/README" , 'w') as outfile:
			outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
	if not os.path.exists(of + "/.internal/parseltongue/README"):	
		with open(of + "/.internal/parseltongue/README" , 'w') as outfile:
			outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")

	# Split positions for filtering, if applicable
	if(p.keep_positions != "none" and os.path.exists(p.keep_positions)):
		click.echo(gettime() + "Found file to filter for specific positions: " + p.keep_positions)
		split_R = script_dir + "/R/01_splitSNPs.R" 
		r_callSplit1 = " ".join(["Rscript", split_R, temp_split, p.keep_positions, "keep"])
		os.system(r_callSplit1)
	
	if(p.remove_positions != "none" and os.path.exists(p.remove_positions)):
		click.echo(gettime() + "Found file to filter against specific positions: " + p.remove_positions)
		split_R = script_dir + "/R/01_splitSNPs.R" #
		r_callSplit2 = " ".join(["Rscript", split_R, temp_split, p.remove_positions, "remove"])
		os.system(r_callSplit2)

	# Do this from above 
	chrs_go = bam_chrs
	
	click.echo(gettime() + "Parsing reads from these chromosomes: ")
	click.echo(chrs_go)

	# Write them to a file
	with open(of + "/.internal/chrs.txt", 'w') as t:
		for item in chrs_go:
			t.write("%s\n" % item)
	
	# Split bams by chromosome and sort by name
	line1 = 'python ' +script_dir+'/python/02_splitBam.py --input '+p.bamfile + " --ncores " + str(ncores)
	line2 =   ' --chrfile ' + of + "/.internal/chrs.txt" + ' --out ' + of + " --barcode-tag " +  p.barcode_tag
	filt_split_cmd = line1 + line2
	
	os.system(filt_split_cmd)
	
	# Let Snakemake process the chromosome files
	y_s = of + "/.internal/parseltongue/hisnpper.object.yaml"
	with open(y_s, 'w') as yaml_file:
		yaml.dump(dict(p), yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
	cp_call = "cp " + y_s +  " " + logs + "/" + p.name + ".parameters.txt"
	os.system(cp_call)
	
	# Setup snakemake logs
	snake_stats = logs + "/" + p.name + ".snakemake_edits.stats"
	snake_log = logs + "/" + p.name + ".snakemake_edits.log"

	exit("exiting before snakemake")
	
	snake_log_out = ""
	if not snake_stdout:
		snake_log_out = ' &>' + snake_log 
		
	snakecmd_chr = 'snakemake --snakefile '+script_dir+'/snake/Snakefile.edits.txt --cores '+str(ncores)+' --config cfp="' + y_s + '" --stats '+snake_stats+snake_log_out
	os.system(snakecmd_chr)

	if keep_temp_files:
		click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.")
	else:
		byefolder = of
		shutil.rmtree(byefolder + "/.internal")
		shutil.rmtree(byefolder + "/temp")
	
		click.echo(gettime() + "Intermediate files successfully removed.")

	click.echo(gettime() + "Complete.")

