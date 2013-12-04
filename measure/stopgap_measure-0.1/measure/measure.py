#!/usr/bin/env python
'''
Methylation estimation of bisulphite-treated DNA.

Usage: measure.py [options] -b <bamfile.bam> -f <fastafile>

Arguments:
 -b, --bam       a bam file of interest. Must end in '.bam'
 -f, --fasta     the accompanying indexed fasta reference file

Optional arguments:
 -m=<mask>       where <mask> is 'cg' or 'c', to only consider
 --mask=<mask>     CpG or CpN sites. 'cg' by default. 
 -s, --sense     only save results for the sense strand reads, 
                   both strands by default.
 -a, --antisense only save results for the antisense strand reads
 -x, --excel     save output in Microsoft Excel format
 -h, --help      show this help message and exit
 -v, --version   prints the version and a disclaimer 

Testing:
 -t, --test      runs self-tests (not yet implemented)
 
Description:
A samfile walker which uses csamtools pileup to add up
(for each unmasked nucleotide) the number of nucleotides, INDELS
and read starts at each position. In addition to presenting this
information, for each position a methylation rate is formed by
dividing C / (C + T).

The results can be exported to an Excel spreadsheet or csv file

To generate the alignments any aligner that generates a legal
BAM file may be used. To index the FASTA reference file,
use 'samtools faidx'.

The output will be saved in the current working directory with
the same name as the BAM file except with 
a .csv or .xlsx extension. 
'''

import sys
import getopt
import re
import pysam
import pandas
import numpy as np
import os
import ctypes
from distutils.version import StrictVersion

#import multiprocessing

version = '0.1'

# http://doeidoei.wordpress.com/2009/03/22/python-tip-3-checking-available-ram-with-python/
class MemoryCheck():
    ''''Checks memory of a given system'''

    def __init__(self):

        if os.name == "posix":
            self.value = self.linux_ram()
        elif os.name == "nt":
            self.value = self.windows_ram()
        else:
            print "I only work with Win or Linux :P"

    def windows_ram(self):
        ''''Uses Windows API to check RAM in this OS'''
        kernel32 = ctypes.windll.kernel32
        c_ulong = ctypes.c_ulong
        class MEMORYSTATUS(ctypes.Structure):
            _fields_ = [
                ("dwLength", c_ulong),
                ("dwMemoryLoad", c_ulong),
                ("dwTotalPhys", c_ulong),
                ("dwAvailPhys", c_ulong),
                ("dwTotalPageFile", c_ulong),
                ("dwAvailPageFile", c_ulong),
                ("dwTotalVirtual", c_ulong),
                ("dwAvailVirtual", c_ulong)
            ]
        memoryStatus = MEMORYSTATUS()
        memoryStatus.dwLength = ctypes.sizeof(MEMORYSTATUS)
        kernel32.GlobalMemoryStatus(ctypes.byref(memoryStatus))

        return int(memoryStatus.dwTotalPhys/1024**2)

    def linux_ram(self):
        """Returns the RAM of a linux system"""
        totalMemory = os.popen("free -m").readlines()[1].split()[1]
        return int(totalMemory)




class ReadMatrix():
    '''
    A samfile walker for methylation estimation of bisulphite-treated DNA.
    
    
    
    class which takes a reference within a samfile
    and tabulates the number of nucleotides, INDELs
    and \'methylation\' rate for each unmasked nucleotide.

    Constructor requires a samfile, a reference
    and optional filters, fasta file and mask.

    '''
    def __init__(self, samfile, reference, start=None, end=None, subset=None, fastafile=None, mask=None, auto=True):
        '''
       Instances are constructed with the below default arguments.

       class ReadMatrix:
            start = None
            end = None
            subset = None
            fastafile = None
            mask = None
            auto = True
        '''
        if isinstance(samfile, pysam.csamtools.Samfile):
            self.sam = samfile
        else:
            self.sam = pysam.Samfile(samfile, 'rb' )
        
        if reference in self.sam.references:
            pass
        else:
            print "Reference not in samfile"
            print samfile.references
            raise ValueError
        
        self.reference = reference
        self.start = start
        self.end = end
        self.subset = subset
        self.length = self.sam.lengths[self.sam.references.index(self.reference)]
        self.fasta = None
        self.mask = None
        
        if fastafile is not None:
            self.load_fasta(fastafile)
        
        if mask is not None:
            if fastafile is None:
                print "You need a reference fasta file to mask the sequence"
                raise ValueError
            if mask != 'c' and mask != 'cg':
                print "Only 'c' and 'cg' are valid masks"
                raise ValueError
            if mask == 'cg':
                self.cg_mask()
            if mask == 'c':
                self.c_mask()
        
        if auto is True:
            self.make_pileup()
    
    def all_counts(self):
        if self.mask is None:
            both = pandas.DataFrame(self._counts[0]) + pandas.DataFrame(self._counts[1])
        else:
            both = pandas.DataFrame(self._counts[0], index = self.mask) + pandas.DataFrame(self._counts[1], index = self.mask)
        #both['Meth_rate'] = self._counts[1]['C'] / (self._counts[1]['C'] + self._counts[1]['T'])
        both['meth_rate'] = self.methylation_rate(both)
        both = pandas.concat([self._ref_bases.ix[:,['reference','seq']], both], axis = 1)
        return both

    def sense_counts(self):
        if self.mask is None:
            sense = pandas.DataFrame(self._counts[0])
        else:
            sense = pandas.DataFrame(self._counts[0], index = self.mask)
        #sense['Meth_rate'] = self._counts[0]['C'] * 1.0 / (self._counts[0]['C'] + self._counts[0]['T'])
        sense['meth_rate'] = self.methylation_rate(sense)
        sense = pandas.concat([self._ref_bases.ix[:,['reference','seq']], sense], axis = 1)
        return sense

    def antisense_counts(self):
        if self.mask is None:
            antis = pandas.DataFrame(self._counts[1])    
        else:
            antis = pandas.DataFrame(self._counts[1], index = self.mask)
        #antis['Meth_rate'] = self._counts[1]['C'] * 1.0 / (self._counts[1]['C'] + self._counts[1]['T'])
        antis['meth_rate'] = self.methylation_rate(antis)
        antis = pandas.concat([self._ref_bases.ix[:,['reference','seq']], antis], axis = 1)
        return antis
    
    def init_counts_table(self):
        
        col_names = ['n', 'A', 'C', 'G', 'T', 'N', 'Del', 'Ins', 'Starts', 'Ends']
        col_dtypes = ['<u4'] * len(col_names)
        col_dict = {'names': col_names, 'formats': col_dtypes }
        
        if self.mask is None:
            self._counts = np.zeros([2, self.length], dtype=col_dict)
        if self.mask is not None:
            self._counts = np.zeros([2, len(self.mask)], dtype=col_dict)
    
    def load_fasta(self, fastafile):
        if isinstance(fastafile, pysam.csamtools.Fastafile):
            self.fasta = fastafile
        else:
            self.fasta = pysam.Fastafile(fastafile)
    
    def regex_mask(self, pattern = ".", select = 1):
        
        seq = self.fasta.fetch(reference=self.reference)
        locations = []
        for m in re.finditer(pattern.lower(), seq.lower()):
            for i in range(m.start(), m.start() + select):
                locations.append(i) 
        self.mask = locations
    
    def c_mask(self):
        '''Mask all bases in the reference except for c or y'''
        self.regex_mask(pattern = 'c|y')
    
    def cg_mask(self):
        '''Mask all bases in the reference except for c or y in the context cg or yg'''
        self.regex_mask(pattern = 'cg|yg')
    
    def remove_mask(self):
        '''Remove an existing reference nucleotide mask'''
        self.mask = None

    def make_pileup(self, max_depth = float('inf'), *args, **kwargs):
        # Make an iterator across the particular reference positions
        self.init_counts_table()
        pileup_iter = self.sam.pileup(reference=self.reference, start=self.start, end=self.end, *args, **kwargs)
        #pileup_iter.addReference()
        # Traverse the iterator by position
        mask_index_start = 0
        for col_region in pileup_iter:
            i = col_region.pos
                      
            if self.mask is not None:
                try:
                    i = self.mask[mask_index_start:].index(i)
                    i =  mask_index_start + i
                    mask_index_start = i
                except ValueError:
                    continue
            
            for pileupread in col_region.pileups:
                '''
                A read flagged as is_del in the pysam documentation is:
                '1 iff the base on the padded read is a deletion'
                So these bases shouldn't be included in the tally as they're
                insertions to the reference.
                '''
                strand = pileupread.alignment.is_reverse
                if pileupread.is_del == 1:
                    self._counts[strand]['Del'][i] += 1
                else:
                    self._counts[strand][pileupread.alignment.seq[pileupread.qpos]][i] += 1
                #indel length; 0 for no indel, positive for ins and negative for del
                if pileupread.indel > 0: #insertion
                    #self._counts['Ins'][i:i + pileupread.indel - 1] += 1
                    self._counts[strand]['Ins'][i] += pileupread.indel
                if pileupread.is_head > 0:
                    self._counts[strand]['Starts'][i] += 1
                if pileupread.is_tail > 0:
                    self._counts[strand]['Ends'][i] += 1

            #for z in [0,1]:
                #self._counts[j]['pos'][i] = col_region.pos
                # Make n consistent with IGV
            #    self._counts[z]['n'][i] =  self._counts[z]['A'][i] + self._counts[z]['C'][i] + self._counts[z]['G'][i] + self._counts[z]['T'][i] + self._counts[z]['N'][i]
        self._counts['n'] =  self._counts['A'] + self._counts['C'] + self._counts['G'] + self._counts['T'] + self._counts['N']

        bases = list(self.fasta.fetch(reference=self.reference))
        if self.mask is None:
            self._ref_bases = pandas.DataFrame({ 'seq': np.array(bases) })
        else:
            self._ref_bases = pandas.DataFrame({ 'seq': np.array(bases)[self.mask] }, index=self.mask)
        self._ref_bases['reference'] = self.reference
        
    def methylation_rate(self, df):
        '''Returns the methylation rate by dividing C / (C + T)'''
        return df['C'] * 1.0 / (df['C'] + df['T'])
    
    def close(self):
        '''Close the bamfile (and fasta file if required)'''
        self.sam.close()
        if self.fasta is not None:
            self.fasta.close()
        

def usage():
    '''Print usage help'''
    print "USAGE:    measure.py  [-hvxsa] -b <bamfile.bam> -f <fastafile> [-m=<mask>]"
    print ""
    print "Walk through a BAM file and generate a methylation report"
    print "For help type: measure.py -h"

def disclaimer():
    '''Print the CSIRO license disclaimer'''
    print 'measure ' + version
    print 'Copyright (C) 2012 CSIRO'
    print 'License CSIRO Open Source Software License, a modification of the GNU GPL version 2 license.'
    print 'This license provides warranties and limitations provisions in addition to those of the GPL.'
    print 'CSIRO license: <http://www.ict.csiro.au/downloads.php>'
    print ''
    print 'Written by Jason Ross'


def main():
    
    bam_file = None
    fasta_file = None
    excel = False
    mask = 'cg'
    strands = 'all'
    results = []

    # parse command line options
    try:
        opts, remaining = getopt.getopt(sys.argv[1:], 'b:f:m:x:sahtv', ['bam=', 'fasta=', 'mask=', 'excel', 'sense', 'antisense', 'help', 'test', 'version'])
    except getopt.error, msg:
        print msg
        usage()
        disclaimer()
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ('-b', '--bam'):
            bam_file = arg
        elif opt in ('-f', '--fasta'):
            fasta_file = arg
        elif opt in ('-m', '--mask'):
            mask = arg
        elif opt in ('-x', '--excel'):
            if StrictVersion(pandas.__version__) >= StrictVersion('0.8.0'):
                excel = True
            else:
                print "Pandas version is " + pandas.__version__ + ", require 0.8.0 or greater for Excel output"
        elif opt in ('-s', '--sense'):
            strands = 'sense'
        elif opt in ('-a', '--antisense'):
            # in the case where both --sense and --antisense are given
            if strands == 'sense':
                strands = 'all'
            else:
                strands = 'antisense'
        elif opt in ('-h', '--help'):
            print __doc__
            sys.exit(0)
        elif opt in ('-t', '--test'):
            print "Not yet implemented"
            disclaimer()
            sys.exit(0)
        elif opt == '--version':
            disclaimer()
            sys.exit(0)

    
    # Sanity checks
    if bam_file is None or fasta_file is None:
        usage()
        sys.exit(0)

    if strands not in ['all', 'sense', 'antisense']:
        print "Illegal strands value " + strands
        sys.exit(2)
    
    if mask not in ['cg', 'c' , 'None']:
        print "Illegal mask value " + mask
        sys.exit(2)
    
    if mask is 'None':
        mask = None
    
    bam_name = re.sub('(.*)\.bam', '\g<1>', bam_file)
    bam_name = bam_name + '_' + strands
    
    if excel is True:
        excel_writer = pandas.ExcelWriter(bam_name + '.xlsx')
    
    csam = pysam.Samfile(bam_file, "rb" )
    
    for ref in csam.references:
        m = ReadMatrix(csam, ref, fastafile=fasta_file, mask=mask)
        if strands == 'all':
            ref_counts = m.all_counts()
        elif strands == 'sense':
            ref_counts = m.sense_counts()
        elif strands == 'antisense':
            ref_counts = m.antisense_counts()
        results.append(ref_counts)    
    
    all_results = pandas.concat(results)
    if excel is True:
        all_results.to_excel(excel_writer, header=True, index=True, index_label='pos')
    else:
        all_results.to_csv(bam_name + '.csv', header=True, index=True, index_label='pos')
    
    # Tidy up and close file handles
    m.close()
    if excel is True:
        excel_writer.save()
    

if __name__ == "__main__":
    main()
