# Authors: Dustin Schones and Keji Zhao
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# schonesde@mail.nih.gov)


"""
This module contains the classes to deal with UCSC
known gene files
"""
import re, os, shutil, time, sys
from math import *
from string import *

plus = re.compile('\+');
minus = re.compile('\-');

import BED;

UCSCError = "Error in UCSC class";


class UCSC:
    """
    Class for keeping known gene information.  This might be too much
    information to carry around when the most usefull parts are really
    the chrom, strand, txStart and txEnd.
    """

    def __init__(self, name, chrom, strand, txStart, txEnd, cdsStart,
                 cdsEnd, exonCount, exonStarts, exonEnds):
        self.name = name;
        self.chrom = chrom;
        self.strand = strand;
        self.txStart = txStart;
        self.txEnd = txEnd;
        self.cdsStart = cdsStart;
        self.cdsEnd = cdsEnd;
        self.exonCount = exonCount;
        self.exonStarts = exonStarts;
        self.exonEnds = exonEnds;
        
    def __set__(self, name, chrom, strand, txStart, txEnd, cdsStart,
                cdsEnd, exonCount, exonStarts, exonEnds):
        self.name = name;
        self.chrom = chrom;
        self.strand = strand;
        self.txStart = txStart;
        self.txEnd = txEnd;
        self.cdsStart = cdsStart;
        self.cdsEnd = cdsEnd;
        self.exonCount = exonCount;
        self.exonStarts = exonStarts;
        self.exonEnds = exonEnds;

    def getAll(self):
        outstring = self.name + " " + self.chrom + " " + self.strand + " " + \
                   str(self.txStart) + \
                    " " + str(self.txEnd) + " " + str(self.cdsStart) + " " + \
                    str(self.cdsEnd) + " " + str(self.exonCount) + " " + \
                    str(self.exonStarts) + " " + str(self.exonEnds);
        try:
            return outstring;
        except:
            sys.stderr.write("No UCSC known gene information for %s\n" % self)
            return ''
    

#----------------------------
#----------------------------

class UCSC_lite:
    """
    Class for keeping known gene information.  Only partial (most
    usefull) information stored.
    """
    def __init__(self, name, chrom, strand, txStart, txEnd):
        self.name = name;
        self.chrom = chrom;
        self.strand = strand;
        self.txStart = txStart;
        self.txEnd = txEnd;
    def __set__(self, name, chrom, strand, txStart, txEnd):
        self.name = name;
        self.chrom = chrom;
        self.strand = strand;
        self.txStart = txStart;
        self.txEnd = txEnd;

    def getAll(self):
        outstring = self.name + " " + self.chrom + " " + self.strand + \
                    " " + str(self.txStart) + " " + str(self.txEnd);
        try:
            return outstring;
        except:
            sys.stderr.write("No UCSC known gene information for %s\n" % self)
            return ''


#----------------------------
#----------------------------


class KnownGenes:
    """
    Class to read in UCSC known gene files and get all the info
    """
    
    def __init__(self, file=None):
        """
        Reads in a UCSC known gene file and builds a dictionary
        with chromosomes as keys and lists of bed vals representing
        the genes on each chromosome as values.
        
        """
        self.gene_coords = {}
        infile = open(file);
        for line in infile:
            """ check to make sure not a header line """
            if not re.match("#", line):
                line = line.strip();
                sline = line.split();
                if sline[1] not in self.gene_coords.keys():
                    self.gene_coords[sline[1]] = [];
                coord = UCSC(sline[0], sline[1], sline[2], atoi(sline[3]),
                             atoi(sline[4]), atoi(sline[5]), atoi(sline[6]),
                             atoi(sline[7]), sline[8], sline[9]);
                self.gene_coords[sline[1]].append(coord);


    def getPromoters(self, upstream, downstream):
        """
        Return the promoter coordinates defined by given upstream and
        downstream distances in a bed dictionary

        NOTICE:  if gene is on negative strand, the promoter start
        and end are still sequential
        """
        self.prom_coords = {};
        for chrom in self.gene_coords.keys():
            chrom_coords = self.gene_coords[chrom];
            for g in chrom_coords:
                if plus.match(g.strand):
                    prom_start = g.txStart - upstream;
                    prom_end = g.txStart + downstream;
                elif minus.match(g.strand):
                    prom_start = g.txEnd - downstream;
                    prom_end = g.txEnd + upstream;
                if prom_start > 0:
                    ucsc = UCSC_lite(g.name, chrom, g.strand,
                                     prom_start, prom_end);
                else:
                    ucsc = UCSC_lite(g.name, chrom, g.strand,
                                     0, prom_end);
                if chrom not in self.prom_coords.keys():
                    self.prom_coords[chrom] = [];
                self.prom_coords[chrom].append(ucsc);
        return self.prom_coords;
        

    def keys(self):
        """
        Return a list of the keys - duplicating the function of a dictionary
        """
        return self.gene_coords.keys()

    def getNumGenes(self):
        """
        Return the number of genes
        """
        num = 0;
        for c in self.gene_coords.keys():
            num += len(self.gene_coords[c]);
        return num;


    def __setitem__(self, name, bedcoord):
        """
        Sets a new gene coord
        """
        self.gene_coords[name] = bedcoord


    def __getitem__(self, name):
        """
        Returns a bed_val indexed by its name or None if no such bed_val exists
        """
        if self.gene_coords.has_key(name):
            return self.gene_coords[name]
        else: raise UCSCError
    

    def __del__(self):
        """
        Delete, delete;
        """
        self.gene_coords.clear()

    def __contains__(self, item):
        """
        Returns  mapping iterator
        """
        return self.gene_coords.has_key(item)

    def __iter__(self):
        """
        Returns mapping iterator
        """
        return self.gene_coords.iterkeys()

    def __len__(self):
        """
        Returns number of gene_coords
        """
        return len(self.gene_coords)

    def __delitem__(self, name):
        """
        removes a chrom if its name exists in the dictionary
        -- I guess this could possibly be useful at some point
        
        """
        if self.gene_coords.has_key(name):
            del self.gene_coords[name]
        else: raise bedError



