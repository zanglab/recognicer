#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import BED
import GenomeData;
import associate_tags_with_regions
import separate_by_chrom
import get_total_tag_counts
import utility
import scipy
from scipy import scipy.stats

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--rawchipreadfile", action="store", type="string", dest="chipreadfile", metavar="<file>", help="raw read file from chip in bed format")
	parser.add_option("-b", "--rawcontrolreadfile", action="store", type="string", dest="controlreadfile", metavar="<file>", help="raw read file from control in bed format")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>", help="average size of a fragment after CHIP experiment")
	parser.add_option("-d", "--islandfile", action="store", type="string", dest="islandfile", metavar="<file>", help="island file in BED format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="island read count summary file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		genomesize = sum (GenomeData.species_chrom_lengths[opt.species].values());
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	chip_library_size=get_total_tag_counts.get_total_tag_counts(opt.chipreadfile);
	control_library_size=get_total_tag_counts.get_total_tag_counts(opt.controlreadfile);
	print "chip library size  ", chip_library_size
	print "control library size  ", control_library_size
	
	totalchip = 0;
	totalcontrol = 0;
	
	islands = BED.BED(opt.species, opt.islandfile, "BED3", 0);
	
	# separate by chrom the chip library
	if utility.fileExists(opt.chipreadfile):
		separate_by_chrom.separateByChrom(chroms, opt.chipreadfile, '.bed1');
	else:
		print opt.chipreadfile, " not found";
		sys.exit(1)
	# separate by chrom the control library
	if utility.fileExists(opt.controlreadfile):
		separate_by_chrom.separateByChrom(chroms, opt.controlreadfile, '.bed2');
	else:
		print opt.controlreadfile, " not found";
		sys.exit(1)	
	
	island_chip_readcount = {};
	island_control_readcount = {};
	
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:
				island_list = islands[chrom];
				if utility.is_sorted(island_list) == 0:
					island_list.sort(key=operator.attrgetter('start'));
					
				island_start_list = []
				island_end_list = []
				for item in island_list:
					island_start_list.append(item.start)
					island_end_list.append(item.end)
	
				island_chip_readcount_list=[0]*len(island_list);
				read_file = chrom + ".bed1";
				f = open(read_file,'r')
				for line in f:
					if not re.match("#", line):
						line = line.strip()
						sline = line.split()
						position = associate_tags_with_regions.tag_position(sline, opt.fragment_size)
						index =associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position);
						if index >= 0:
							island_chip_readcount_list[index] += 1;
							totalchip += 1;
				f.close();
				island_chip_readcount[chrom] = island_chip_readcount_list;
							
				island_control_readcount_list=[0]*len(island_list);
				read_file = chrom + ".bed2";
				f = open(read_file,'r')
				for line in f:
					if not re.match("#", line):
						line = line.strip()
						sline = line.split()
						position = associate_tags_with_regions.tag_position(sline, opt.fragment_size)
						index = associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position);
						if index >= 0:
							island_control_readcount_list[index] += 1;
							totalcontrol += 1;
				f.close();
							
				island_control_readcount[chrom] = island_control_readcount_list;			
						
	chip_background_read = chip_library_size - totalchip;
	control_background_read = control_library_size - totalcontrol;
	#scaling_factor = chip_background_read*1.0/control_background_read;
	scaling_factor = chip_library_size*1.0/control_library_size;

	print "Total number of chip reads on islands is: ", totalchip; 
	print "Total number of control reads on islands is: ", totalcontrol; 

	print "chip_background_read   ", chip_background_read
	print "control_background_read   ", control_background_read

	out = open(opt.out_file, 'w');
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:					
				island_list = islands[chrom];			
				for index in xrange(len(island_list)):
					if (island_control_readcount[chrom])[index] > 0:
						average = (island_control_readcount[chrom])[index] * scaling_factor;
					else:
						length = item.end-item.start;		average = length * control_library_size *1.0/genomesize;			
						average = min(0.25, average)* scaling_factor;
					observation = (island_chip_readcount[chrom])[index];
					if observation > average:
						pvalue =	 scipy.stats.poisson.sf((island_chip_readcount[chrom])[index], average)[()]; 
					else:
						pvalue = 1;
					
					item = island_list[index];
					outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + str(observation) + "\t" + str((island_control_readcount[chrom])[index]) + "\t" + str(pvalue) + "\n";	
					out.write(outline);		
	out.close();
	
	
	separate_by_chrom.cleanup(chroms, '.bed1');
	separate_by_chrom.cleanup(chroms, '.bed2');
	
	

if __name__ == "__main__":
	main(sys.argv)
