import click
import os
import os.path
import sys
import shutil
import random
import string
import itertools
import time
import platform
from ruamel import yaml
from pkg_resources import get_distribution

from .frombapHelp import *

class hisnpperProject():
	def __init__(self, script_dir, mode,
			bamfile, snps, fasta, output,
			name, barcode_tag, haplotype_tag, min_aq,
			ncores, keep_temp_files):
			
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.hisnpper_version = get_distribution('hisnpper').version
		self.script_dir = script_dir
		self.mode = mode
		self.bamfile = bamfile
		self.snps = snps
		self.fasta = fasta
		self.output = output
		self.name = name			
		self.barcode_tag = barcode_tag
		self.haplotype_tag = haplotype_tag
		self.min_aq = str(min_aq)
		self.ncores = ncores
				
		if(name == "default"):
			filename, file_extension = os.path.splitext(self.bamfile)
			self.name = os.path.basename(filename)
		
	#--------------------------------------------------------------------------------
	# Define a method to dump the object as a .yaml/dictionary for use in other files
	#--------------------------------------------------------------------------------
	def __iter__(self):
		yield 'hisnpper_version', self.hisnpper_version
		yield 'script_dir', self.script_dir
		yield 'mode', self.mode
		yield 'bamfile', self.bamfile
		yield 'snps', self.snps
		yield 'fasta', self.fasta
		yield 'output', self.output
		yield 'name', self.name			
		yield 'barcode_tag', self.barcode_tag
		yield 'haplotype_tag', self.haplotype_tag
		yield 'min_aq', self.min_aq
		yield 'ncores', self.ncores

class hisnpperEditProject():
	def __init__(self, script_dir,
			bamfile, fasta, name, output,
			edit, min_aq, barcode_tag, keep_positions, remove_positions,
			ncores, keep_temp_files):
		
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.hisnpper_version = get_distribution('hisnpper').version
		self.script_dir = script_dir
		self.bamfile = bamfile
		self.fasta = fasta
		self.name = name	
		self.output = output
		self.edit = edit

		self.min_aq = str(min_aq)
		self.barcode_tag = barcode_tag
		self.keep_positions = keep_positions
		self.remove_positions = remove_positions
		self.ncores = ncores
				
		if(name == "default"):
			filename, file_extension = os.path.splitext(self.bamfile)
			self.name = os.path.basename(filename)
		
	#--------------------------------------------------------------------------------
	# Define a method to dump the object as a .yaml/dictionary for use in other files
	#--------------------------------------------------------------------------------
	def __iter__(self):
		yield 'hisnpper_version', self.hisnpper_version
		yield 'script_dir', self.script_dir
		yield 'bamfile', self.bamfile
		yield 'fasta', self.fasta
		yield 'name', self.name
		yield 'output', self.output
		yield 'edit', self.edit

		yield 'min_aq', self.min_aq	
		yield 'barcode_tag', self.barcode_tag		
		yield 'keep_positions', self.keep_positions
		yield 'remove_positions', self.remove_positions

		yield 'ncores', self.ncores

