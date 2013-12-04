#!/usr/bin/env python
'''
Created on 17/06/2013

@author: ros259
'''

version = 0.01
debug = True

import multiprocessing as mp
import pysam


#from scipy.stats import itemfreq
import itertools
import warnings
import sys
import re
from realign import *
import traceback

#import pdb

#import copy_reg
#import types


'''
def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    if func_name.startswith('__') and not func_name.endswith('__'):
        cls_name = cls.__name__.lstrip('_')
        func_name = '_' + cls_name + func_name
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.__mro__:
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)


copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


def _map_unwrap(arg, **kwarg):
 
    try:
        realigned = RealignReads._RealignReadTest(*arg, **kwarg)
    except:
        raise "Multiprocessor exception"
    return realigned
    return RealignReads._RealignReadTest(*arg, **kwarg)
'''

class RealignReads(object):
    '''
    ref: dict
    read: SeqRecord
    '''
    def __init__(
                    self, sam_in=None, ref=None, sam_out=None, matrix=None,
                    gap_open=-16, gap_extend=-4, gap_char='-',
                    reverse_align=True, only_gapped=False, compute_scores=True,
                    binary_mode=True,
                    softclip_limit=5,
                    cpus=mp.cpu_count()
                ):
        if sam_in is not None and not isinstance(sam_in, pysam.Samfile):
            raise AttributeError('Need pysam samfile object')
        if ref is not None and not isinstance(ref, pysam.Fastafile):
            raise AttributeError('Need pysam samfile object')

        self.sam_in = sam_in
        self.ref = ref
        self.sam_out = sam_out
        if matrix is None:
            self.matrix = 'bismat.txt'
        else:
            self.matrix = matrix
        self.only_gapped = only_gapped
        self.binary_mode = binary_mode
        self.softclip_limit = softclip_limit
        self.cpus = cpus
        if sam_in is None:
            self.refnames = None
        else:
            self.refnames = sam_in.references

        self.matrix_scores = ScoringMatrix(filename=self.matrix)
        self.realign = Realign(
            matrix=self.matrix, gap_open=gap_open,
            gap_extend=gap_extend, gap_char=gap_char,
            reverse_align=reverse_align, compute_score=compute_scores
        )

    def _find_reference_coords(self, read):

        is_read_start = True
        ref_slice = {'start': read.pos, 'end': read.pos + len(read.seq)}   # is 0-based

        for c in read.cigar:
            '''
            Need to adjust reference start and stop for softclips and ref
            insertions.
            Softclips may only have hardclip operations between them and
            the ends of the CIGAR.
            '''
            if c[0] == 1:    # insertion to the ref
                ref_slice['end'] -= c[1]
            if c[0] == 2:    # deletion from the ref
                ref_slice['end'] += c[1]
            elif c[0] == 3:  # skipped from the ref is not defined here
                warnings.warn('Reads with N in the CIGAR is not defined')
            elif c[0] == 4:  # softclip
                if is_read_start is True:
                    ref_slice['start'] -= c[1]
                else:
                    ref_slice['end'] += c[1]
            if is_read_start is True and c[0] != 5:
                # Special case for hardclips
                is_read_start = False

        if ref_slice['start'] < 0:
            ref_slice['start'] = 0

        return(ref_slice)

    def Remap(self, count=None, verbose=False, debug=None):
        #scores = {}
        if self.sam_in is None:
            return
        counter = 0
        has_score = False

        if self.sam_out is None:
            self._out = sys.stdout
        else:
            write_mode = 'wb'
            if self.binary_mode is False:
                write_mode = 'wh'
            self._out = pysam.Samfile(
                self.sam_out,
                mode=write_mode,
                referencenames=self.sam_in.references,
                referencelengths=self.sam_in.lengths,
                header=self._MakeHeader(self.sam_in.header)
            )

        for counter, read in enumerate(self.sam_in.fetch()):
            'Optional setting of count, to only realign count reads'
            if count is not None and counter > count:
                break
            if read.is_unmapped is True:
                self._out.write(read)
                continue
            if counter == 0:
                # Check if an alignment score is already present
                # If it is then record this in the has_score flag
                tags = read.tags
                for i in range(0, len(tags)):
                    if tags[i][0] == 'AS':
                        has_score = True
                        continue

            ''' If we wish to only realign gapped reads then check each read
            to see if it has a gap. If it does, then write the read to output
            and continue to the next read'''
            has_indel = False
            num_softclips = sum(
                [length for c, length in [read.cigar[0], read.cigar[-1]]
                if c == 4])
            #print('%s -> %s' % (read.cigar, num_softclips))
            for c in read.cigar:
                if c[0] == 1 or c[0] == 2:
                    # read has an indel
                    has_indel = True
                    break

            if self.only_gapped is True and has_indel is False:
                self._out.write(read)
                #if the read is a perfect match then don't realign
                continue

            ref_slice = self._find_reference_coords(read)
            ref = self.ref.fetch(
                reference=self.refnames[read.tid],
                start=ref_slice['start'], end=ref_slice['end']
            )

            '''
            Realignment, scoring, new position and cigar
            '''
            
            '''
            if num_softclips < self.softclip_limit:
                aligntype='global'
            else:
                aligntype='local'
            try:
                realn = self.realign.align(
                    ref=ref,
                    query=read.seq,
                    r_offset=read.pos,
                    aligntype=aligntype,
                    ref_name=self.refnames[read.tid],
                    query_name=read.qname
                )
            except Exception, e:
                print("Failed on read %s" % read.qname)
                sys.stderr.write(repr(e))
                exit()
            '''
            try_align = True
            aligntype = 'global'
            while try_align is True:
                try:
                    realn = self.realign.align(
                        ref=ref,
                        query=read.seq,
                        #r_offset=read.pos,
                        r_offset=ref_slice['start'],
                        aligntype=aligntype,
                        ref_name=self.refnames[read.tid],
                        query_name=read.qname
                    )
                except Exception, e:
                    print("Failed on read %s" % read.qname)
                    sys.stderr.write(repr(e))
                    exit()
                #if aligntype == 'global' and realn.score < 0:
                    #aligntype = 'local'
                else:
                    try_align = False

            if aligntype == 'global' and realn.score < 0:
                # give up
                self._out.write(read)
                continue

            if aligntype == 'local' and realn.score < 0:
                print('Bad alignment for %s' % read.qname)

            if verbose is True:
                    realn.dump()

            if debug is not None:
                for read_string in debug:
                    r = re.search(read_string, read.qname)
                    if r is not None:
                        print('################   DEBUG   ################')
                        realn.dump()
                if aligntype=='local':
                    print('#############  ALIGN COMPARE  ###############')
                    realn.dump()
                    print("Fixed: {0}".format(MakeSamCigar(SoftClipFix(realn))))
                    realn_global = self.realign.align(
                        ref=ref,
                        query=read.seq,
                        r_offset=ref_slice['start'],
                        #r_offset=read.pos,
                        aligntype='global',
                        ref_name=self.refnames[read.tid],
                        query_name=read.qname
                        )
                    realn_global.dump()                        

            if aligntype == 'local':
            # dirty hack for swalign not supporting softclips
                realn.cigar = SoftClipFix(realn)

            if realn.score is not None:
                if has_score is True:
                    old_score_index = None
                    for i, tag in enumerate(read.tags):
                        if tag[0] == 'AS':
                            old_score_index = i
                            break
                    if old_score_index is None:
                        raise ValueError(
                            "Read {0} is missing an alignment score.".format(read.qname))
                    read.tags[old_score_index] = ('AS', realn.score)
                else:
                    read.tags = [('AS', realn.score)] + read.tags

            read.tags = read.tags + [('OC', MakeSamCigar(realn.cigar)),
                                     ('OP', read.pos + 1)]

            read.cigar = realn.cigar
            read.pos = realn.r_offset + realn.r_pos
            self._out.write(read)
            #counter += 1
        self._out.close


    def PrettyPrint(self, read, realn):
        # Under SAM spec the pos is a 1-based leftmost mapping.
        # POS == 0 is for unmapped reads
        strand = ('+', '-')[read.is_reverse]

        print(read.qname + ' aligns to ' + self.refnames[read.tid]) +\
        ', score(' + str(realn.score) +\
        ')   ####################################'
        print('Read  ' + read.seq)
        print('Sub   ' + realn.aln[1])
        print('Qry   ' + realn.aln[0])
        print(strand + '     ' + str(read.pos + 1) + '  ' +
        str(read.cigar) + '  ' + MakeSamCigar(read.cigar))
        print('->    ' + str(realn.pos + 1) + '  ' +
        str(realn.bam_cigar) + '  ' + MakeSamCigar(realn.bam_cigar))
        

    def MakeBam(self, ref_name, aln):
        qual_list = self.sam_quality()
        qual = ''.join(qual_list[aln.begin:aln.end])
        bam = pysam.AlignedRead()
        bam.qname = self.read.name
        bam.seq = self.read.seq.tostring()
        bam.flag = 0
        #bam.rname = ref_name
        bam.pos = aln.begin + 1
        bam.mapq = 255
        bam.cigar = aln.bam_cigar()
        bam.qual = qual
        bam.tags = (("NM", 1), ("RG", "L1"))
        return bam

    def _MakeHeader(self, header):

        if not isinstance(header, dict):
            raise TypeError("Header needs to be a dict")
            print(header)
            header = {}
        new_pg = {'CL': ' '.join([str(arg) for arg in sys.argv]),
                  'ID': 'stopgap',
                  'PN': 'stopgap',
                  'VN': version}
        #TODO: pysam seems to remove PP tags?
        if 'PG' in header:
            pg_list = header['PG']
            new_pg['PP'] = pg_list[len(pg_list) - 1]['ID']
        else:
            pg_list = []
        pg_list.append(new_pg)
        header['PG'] = pg_list
        return header


        
        #print self.cpus
#    def __call__(self, read):
#        return self._RealignRead(read)
    '''
    def __getstate__(self):
        pass
        return()
    
    def __setstate__(self, state):
        pass
    '''
    
    def _RealignReadList(self, readlist):
        return [self._RealignReadTest(read) for read in readlist]

    def _RealignReadTest(self, read):
        print "Hello"
        return

    def RemapReads(self, count=None):

        counter = 0
        write_mode = 'wb'
        if self.binary_mode is False:
            write_mode = 'wh'
        self._out = pysam.Samfile(self.sam_out,
                                  mode=write_mode,
                                  referencenames=self.sam_in.references,
                                  referencelengths=self.sam_in.lengths,
                                  header=self._MakeHeader(self.sam_in.header)
                                  )
        
        read_iter = itertools.islice(self.sam_in.fetch(), count)        
        p = mp.Pool(self.cpus)

        chunk_size = 10
        while True:        
            aligned_reads = list(itertools.islice(read_iter, chunk_size))
            #aligned_reads = itertools.islice(read_iter, chunk_size)
            # Check if there are any reads left
            print len(aligned_reads)
            if len(aligned_reads) == 0:
                break

            remapped_reads = p.map(_map_unwrap, zip([self]*len(aligned_reads), aligned_reads))
            for read in remapped_reads:
                self._out.write(read)
            
        self._out.close
        return

