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
from .frombapHelp import *

class aaProject():
	def __init__(self, script_dir, input, snp_table, output,
		name, haplotype_tag,
		ncores, keep_temp_files):
		
				
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.script_dir = script_dir
		self.output = output
		self.name = name			
		self.bamfile = input
		self.ncores = ncores
		self.haplotype_tag = haplotype_tag		

		
	#--------------------------------------------------------------------------------
	# Define a method to dump the object as a .yaml/dictionary for use in other files
	#--------------------------------------------------------------------------------
	def __iter__(self):
		
		yield 'script_dir', self.script_dir
		yield 'output', self.output
		yield 'bamfile', self.bamfile
		yield 'name', self.name
		yield 'ncores', self.ncores
		yield 'haplotype_tag', self.haplotype_tag
