#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


import BED
import UCSC
import GenomeData;

grep = "/bin/grep";
cat = "/bin/cat";

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

def fileExists(f):
	try:
		file = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists

def is_bed_sorted(list):
        """
        Check if sorted in ascending order.
        input is a list of BED with chrom, start, end and value.
        output: sorted =1 or 0
        """
	sorted = 1;
	for index in range(0, len(list)-1):
		if list[index].start > list[index+1].start:
			sorted = 0;
	return sorted;
