#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import GenomeData;
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--islandsummary", action="store", type="string", dest="islandsummary", metavar="<file>", help="island summary file")
	parser.add_option("-p", "--pvalue", action="store", type="float", dest="pvalue", metavar="<float>", help="pvalue")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="significant islands under pvalue")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
		
	inputfile = open(opt.islandsummary,'r');
	outfile = open(opt.out_file, 'w');
	total =0;
	for line in inputfile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if atof(sline[5]) <= opt.pvalue:
				total += atoi(sline[3]);
				outfile.write('\t'.join(sline)+'\n');
	inputfile.close()
	outfile.close()
	
	print "Total number of tags on islands with p value ", opt.pvalue, " is ",  total;
	
if __name__ == "__main__":
	main(sys.argv)