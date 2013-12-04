# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 09:50:29 2013

@author: ros259
"""

import sys
import os
import pysam
import cProfile, pstats, io
import pwd
import time

u = pwd.getpwuid(os.getuid())
u.pw_name
base_path = '/home/' + u.pw_name + '/Dropbox/stopgap-measure/'

#cd /home/jason/Dropbox/stopgap-measure
import stopgap

if __name__ == '__main__':
    pr = cProfile.Profile()
    pr.enable()
    bam_file = os.path.join(base_path, 'normal_sorted.bam')    
    ref_fasta_file = os.path.join(base_path, '454BisAmpliconsShrimp.fa')
    samfile = pysam.Samfile(bam_file, "rb")
    reference = pysam.Fastafile(ref_fasta_file)

    test_file = 'normal_stopgap'

    t = stopgap.RealignReads(sam_in=samfile,
                             ref=reference,
                             sam_out=os.path.join(base_path, test_file + '.sam'),
                             matrix=os.path.join(base_path, 'bismat.txt'),
                             gap_open=-10,
                             gap_extend=-6,
                             reverse_align=True,
                             only_gapped=False,
                             compute_scores=True,
                             binary_mode=False,
                             cpus=2)

    #old_sys_stdout = sys.stdout
    #sys.stdout = open('test.txt', 'w')
    t.Remap(verbose=False)
    #t.Remap(2000, verbose=False, debug=['BEFQC', 'CBX7A', 'BB4X9', '1COIY', 'DPNZ4'])
    #GRASPF3R3 is really bad for misalignments due to long softclips
    #t.RemapReadsSingle(verbose=False, debug=['CBX7A', 'BB4X9', '1COIY', 'DPNZ4',
    #'A1QVL', 'BVQ14', 'DJHAG', 'B979S', 'AUWH4', 'B11P8', 'ATXRQ', 'BTLUJ',
    #'A0N2R', 'EGEUA'])
    #sys.stdout.close()
    #sys.stdout = old_sys_stdout
    pr.disable()
    
    #s = io.StringIO()
    ps = pstats.Stats(pr, stream=sys.stdout)
    #ps = pstats.Stats(pr, stream=s)
    ps.sort_stats('time')
    ps.print_stats()
    #print( os.getcwd())
    #time.sleep(3)
    #os.system("samtools view -Sb /home/jason/Dropbox/stopgap-measure/" + test_file + ".sam > /home/jason/Dropbox/stopgap-measure/" + test_file + ".bam")
    #os.system("samtools sort /home/jason/Dropbox/stopgap-measure/" + test_file + ".bam /home/jason/Dropbox/stopgap-measure/" + test_file + "_sorted")
    #os.system("samtools index /home/jason/Dropbox/stopgap-measure/" + test_file + "_sorted.bam")
    #ps.dump_stats('/home/ros259/Dropbox/stopgap-measure/profiler_results.txt')
'''
r = samfile.next()
reference.fetch(reference=samfile.references[r.tid],
                start=r.pos,
                end=r.aend)
'''         
