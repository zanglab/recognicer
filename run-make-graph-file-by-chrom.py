#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

## get BED module
import BED
import GenomeData


"""
For now, we'll use fragment sizes of 150 bp for histone modifications
and 300 bp for enzymes

and window sizes of 200 bp for histone modifications and 400 bp for
enzymes
"""


Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";


#SratDir = os.path.expandvars('$SRAT')
SratDir ='/nv/vol190/zanglab/yw3t/CG9/Modules/'
make_graph_file = os.path.join(SratDir, "make-graph-file.py");


def separateByChrom(chroms, file):
    for chrom in chroms:
        match = chrom + "[[:space:]]";
        tmpFile = chrom + ".bed";
        try:
            if os.system('%s %s %s > %s' %
                         (grep, match, file, tmpFile)): raise
        except: sys.stderr.write("grep failed\n");


def makeGraphFile(chroms, window, fragment_size):
    for chrom in chroms:
        bed_file = chrom + ".bed";
        graph_file = chrom + ".graph";
        try:
            if os.system('%s -f %s -c %s -w %s -i %s -o %s'%
                         (make_graph_file, bed_file, chrom,
                          window, fragment_size, graph_file)): raise
        except: sys.stderr.write("make_graph_file failed\n")


def combineAllGraphFiles(chroms, final_out):
    """
    Combine the seperately processed chromosomes, return the output file name
    """
    outfile = open(final_out,'w');
    outfile.close();
    
    for chrom in chroms:
        graph_file = chrom + ".graph";
        try:
            if os.system('%s %s >> %s' %
                         (cat, graph_file, final_out)): raise
        except: sys.stderr.write("cat failed\n")
    return final_out
		

def cleanup(chroms):
    for chrom in chroms:
        bed_file = chrom + ".bed";
        graph_file = chrom + ".graph";
        try:
            if os.remove('%s' %
                         (bed_file)): raise
        except: sys.stderr.write("clean up failed\n");
        try:
            if os.remove('%s' %
                         (graph_file)): raise
        except: sys.stderr.write("clean up failed\n");


    
def main(argv):
    """
    Note the window_size and the fragment_size are both input as strings, as they are used in
    a shell script in makeGraphFile.
    """
    parser = OptionParser()
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="mm8,hg18,dm2,etc", metavar="<str>")
    parser.add_option("-b", "--bed_file", action="store", type="string",
                      dest="bedfile", help="bed file to make graph file of",
                      metavar="<file>")
    parser.add_option("-w", "--window_size", action="store", type="string",
                      dest="window_size", help="window size", metavar="<int>")
    parser.add_option("-i", "--fragment_size", action="store", type="string",
                      dest="fragment_size",
                      help="size of fragments after CHIP experiment",
                      metavar="<int>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output bed summary file name",
                      metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)

    if opt.species in GenomeData.species_chroms.keys():
        chroms = GenomeData.species_chroms[opt.species];
        separateByChrom(chroms, opt.bedfile);
        makeGraphFile(chroms, opt.window_size, opt.fragment_size);
        final_output_file = opt.outfile;
        final_output_file = combineAllGraphFiles(chroms, final_output_file);
        #in_filename = (opt.bedfile).split('/');
        #in_filename = (in_filename[-1]).split('.');
        #out_filename = in_filename[-2] + "_summary.graph";
        cleanup(chroms);
    else:
        print opt.species + " is not in the species list ";
	
    

if __name__ == "__main__":
	main(sys.argv)


        
