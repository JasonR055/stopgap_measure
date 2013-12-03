# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 18:02:24 2013

@author: ros259
"""
import unittest
import nwalign as nw
import swalign as sw


op_alpha = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7,'X': 8}
op_number = {v: k for k, v in op_alpha.items()}


def MakeSamCigar(bam_cigar):
    sam_cigar = ''
    for op, op_length in bam_cigar:
        #sam_cigar += str(op_length) + 'MIDNSHP=X'[op]
        sam_cigar += str(op_length) + op_number[op]
    return(sam_cigar)


def SoftClipFix(swaln):
    cigar_read_len = 0
    read_length = len(swaln.orig_query)
    new_cigar = []
    for c in swaln.cigar:
        if c[1] in ['M', 'I', 'S', '=', 'X']:
            cigar_read_len += c[0]
        new_cigar.append((op_alpha[c[1]], c[0]))
    if swaln.q_pos > 0:
        new_cigar = [(4, swaln.q_pos)] + new_cigar
        cigar_read_len += swaln.q_pos
    if cigar_read_len < read_length:
        new_cigar = new_cigar + [(4, read_length - cigar_read_len)]
    return new_cigar


class ScoringMatrix(object):
    '''
    Read an alignment scoring matrix from a file

    The matrix should be space-delimited. For example:

    # Bisulphite tuned matrix
       A  R  N  C  G  T  Y  *
    A  2  2  0 -5 -5 -5 -5 -8
    R  2  2  0 -5  2 -5 -5 -8
    N  0  0  0  0  0  0  0 -8
    C -5 -5  0  2 -5 -5  2 -8
    G -5  2  0 -5  2 -5 -5 -8
    T -5 -5  0 -5 -5  2  2 -8
    Y -5 -5  0  2 -5  2  2 -8
    * -8 -8 -8 -8 -8 -8 -8  1

    The number of matrix rows should equal the number of columns.
    The rows and columns are not required to be ordered.
    Empty lines,or those starting with # are ignored.
    '''

    def __init__(self, filename=None):
        assert filename
        matrix_handle = open(filename)
        self._scores = {}
        self._bases = None
        self._base_count = 0
        row_count = 0
        for l in matrix_handle:
            if l.startswith('#'):
                continue
            if l.startswith('\n'):
                continue
            row = l.split()
            if self._bases is None:
                self._bases = row
                self._base_count = len(self._bases)
                continue
            row_base = row.pop(0)
            if len(row) != self._base_count:
                raise TypeError("Matrix row elements different to alphabet \
                length in cols(" + str(self._base_count) + ")")
            for i in range(self._base_count):
                self._scores[(self._bases[i], row_base)] = int(row[i])
            row_count += 1
        if row_count != self._base_count:
            raise TypeError("Number of matrix rows (" +
            str(row_count) + ") not equal to columns (" +
            str(self._base_count) + ").")

    @property
    def scores(self):
        return self._scores

    @property
    def bases(self):
        return self._bases

    @property
    def base_count(self):
        return self._base_count


class Realign(object):
    '''
    Operation to BAM mapping
    Op # Desc
    M  0 Match or mismatch
    I  1 Insertion to the ref
    D  2 Deletion from the ref
    N  3 Skipped from the ref
    S  4 Soft clipping (clipped sequences present in SEQ)
    H  5 Hard clipping (clipped sequences not present in SEQ)
    P  6 Padding (silent deletion from padded reference)
    =  7 sequence match
    X  8 sequence mismatch

    H needs to be the first or last operation
    Sum of lengths of M/I/S/=/X shall equal the length of SEQ
    '''

    def __init__(self, matrix, gap_open, gap_extend, gap_char,
                 reverse_align=False, compute_score=True):

        self.matrix = matrix
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.gap_char = gap_char
        self.reverse_align = reverse_align
        self.compute_score = compute_score
        self._globalaligninstance = self._globalaligninstance()
        self._localaligninstance = self._localaligninstance()

        # Not using the swalign ScoringMatrix class unless we have to as it's slow

    def align(self, ref, query, r_offset=0, aligntype='global', **kwargs):
        assert(aligntype in ['global', 'local'])
        if aligntype is 'global':
            # An instance of Needle
            aln = self._globalaligninstance.align(ref, query, r_offset, **kwargs)
        else:
            # An instance of swalign.Alignment
            '''
            if self.reverse_align is True:
                orig_query = query
                orig_ref = ref
                query = query[::-1]
                ref = ref[::-1]
            
            aln.cigar.reverse()
            aln.orig_ref = orig_ref
            aln.orig_query = orig_query
            '''
            aln = self._localaligninstance.align(ref, query, **kwargs)            
            aln.set_ref_offset(ref=aln.r_name, offset=r_offset, region=aln.r_region)
        return aln


    def _localaligninstance(self):

        localalign = sw.LocalAlignment(scoring_matrix=sw.ScoringMatrix(self.matrix),
                           gap_penalty=self.gap_open,
                           gap_extension_penalty=self.gap_extend,
                           gap_extension_decay=0.0,
                           prefer_gap_runs=True,
                           verbose=False,
                           globalalign=False)

        return localalign

    def _globalaligninstance(self):

        globalalign = Needle(matrix=self.matrix,
                             gap_open=self.gap_open,
                             gap_extend=self.gap_extend,
                             gap_char=self.gap_char,
                             reverse_align=self.reverse_align,
                             compute_score=self.compute_score)

        return globalalign


class Needle(object):
    
    def __init__(self, matrix, gap_open=-16, gap_extend=-4,
                 gap_char='-', reverse_align=False, compute_score=False):

        #assert(isinstance(scoringmatrix, realign.ScoringMatrix))
        
        self._matrix = matrix
        self._scoringmatrix=ScoringMatrix(matrix)
        self.gap_open=gap_open
        self.gap_extend=gap_extend
        self.gap_char = gap_char
        self.reverse_align = reverse_align
        self.compute_score = compute_score

    @property
    def matrix(self):
        return self._matrix

    @matrix.setter
    def matrix(self, filename):
        self._matrix = filename
        self._scoringmatrix = ScoringMatrix(filename)

    @property
    def scoringmatrix(self):
        return self._scoringmatrix

    def align(self, ref, query, r_offset=0, rc=False, **kwargs):

        orig_ref = ref
        orig_query = query
        query = query.upper()  # index 0
        ref = ref.upper()      # index 1
        if kwargs['query_name']:
            query_name = kwargs['query_name']
        else:
            query_name = 'Query'

        if kwargs['ref_name']:
            ref_name = kwargs['ref_name']
        else:
            ref_name = 'Query'

        score=None
        cigar_list = None
        op_state = None
        op_length = 0
        mismatch_run = 0
        first_match_pos = None
        ref_offset = 0
        
        if self.reverse_align is True:
            query = query[::-1]
            ref = ref[::-1]
        try:
            aln = nw.global_align(query, ref,
                                  gap_open=self.gap_open,
                                  gap_extend=self.gap_extend,
                                  matrix=self.matrix)
        except e:
            print('Query: %s Ref: %s, open:%i extend:%i' %
            (query, ref, self.gap_open, self.gap_extend))
            print('Matrix file: %s' % self.matrix)
            print e

        if self.reverse_align is True:
            aln = tuple(a[::-1] for a in aln)
        #TODO: Clean up cigar code
        for i in range(len(aln[0])):
            op = -1
            r = aln[1][i]
            q = aln[0][i]

            if r == self.gap_char:
                op = 1  # Insertion
                #mismatch_run = 0
            elif q == self.gap_char:
                op = 2  # Deletion
            else:
                op = 0
                this_score = 0
                try:
                    this_score = self.scoringmatrix.scores[(r, q)]
                except KeyError:
                    raise KeyError("Can't find " + r + "/" + q +
                        " in matrix")
                if this_score > 0:
                    mismatch_run = 0
                    if(first_match_pos is None):
                       first_match_pos = i
                if this_score <= 0:
                    mismatch_run += 1

            '''
           Start: Special case for the start of the read
            '''
            if cigar_list is None:
                '''
                Don't count extra reference at the start in the cigar (op=2),
                we need to move the pos instead.
                If there is a mismatch (0 or less match score) or an insertion
                at the start of the read then softclip it.
                '''
                if op == 2:
                    continue
                if op == 1:
                    ref_offset += 1
                if op == 1 or mismatch_run > 0:
                    op = 4
                op_state = op
                op_length = 1
                cigar_list = []
                continue
            '''
           End: Special case for the start of the read
            '''
            if op_state == 4 and op == 2:
                continue

            if op_state == 4 and (mismatch_run > 0 or op == 1):
                op = 4

            # Add to cigar now...
            if op_state == op:
                op_length += 1
            else:
                cigar_list.append((op_state, op_length))
                op_state = op
                op_length = 1

        '''
        Add last element to cigar
        Convert an insertion at the end of the ref into a read softclip
        A deletion from the ref at the end means the ref is too long
        Ignore adding this to the cigar string
        '''

        if op_state == 0 and mismatch_run == 0:
            cigar_list.append((op_state, op_length))
        elif op_state == 0 and mismatch_run > 0:
            cigar_list.append((op_state, op_length - mismatch_run))
            cigar_list.append((4, mismatch_run))
        elif op_state == 1:
            cigar_list.append((4, op_length))
        elif op_state == 2:
            if mismatch_run > 0:
                # TODO:
                # pop off mismatching bases before del and convert to softclip
                prior_cig = cigar_list.pop()
                cigar_list.append((prior_cig[0], prior_cig[1] - mismatch_run))
                cigar_list.append((4, mismatch_run))
            else:
                pass

        cigar_len = 0
        for c in cigar_list:
            if c[0] in [0, 1, 4]:
                cigar_len += c[1]
        if cigar_len != len(query):
            raise ValueError("Cigar length (%r) and \
            read length (%r) differ" % (cigar_len, len(query)))

        if self.compute_score is True:
            #TODO: resolve for softclips
            score = nw.score_alignment(
                aln[0], aln[1], gap_open=self.gap_open,
                gap_extend=self.gap_extend, matrix=self.matrix
            )
        alignment = Alignment(query=query, ref=ref,
                              orig_query=orig_query, orig_ref=orig_ref,
                              aln=aln,
                              q_pos=first_match_pos,
                              r_pos=first_match_pos - ref_offset,
                              r_offset=r_offset,
                              cigar=cigar_list,
                              score=score,
                              ref_name=ref_name,
                              query_name=query_name,
                              rc=rc, globalalign=True)
        return alignment


class Alignment(object):
    def __init__(self, query, ref, orig_ref, orig_query, aln, q_pos, r_pos,
                 r_offset, cigar, score, ref_name='', query_name='',
                 rc=False, globalalign=True):
        
        self.query = query
        self.ref = ref
        self.orig_query = orig_query
        self.orig_ref = orig_ref
        self.aln = aln,
        self.q_pos = q_pos
        self.r_pos = r_pos
        self.r_offset = r_offset
        self.cigar = cigar
        self.score = score
        self.r_name = ref_name
        self.q_name = query_name
        self.rc = rc
        self.globalalign = globalalign

    def dump(self):
        # Under SAM spec the pos is a 1-based leftmost mapping.
        # POS == 0 is for unmapped reads
        aln_string = ''
        for c in self.cigar:
            if c[0] == 0:
                aln_string += '|' * c[1]
            elif c[0] == 4:
                aln_string += '.' * c[1]
            else:
                aln_string += ' ' * c[1]

        #print(self.aln)
        #print(len(self.aln))
        #print(len(self.aln[0]))
        #print(self.aln[0])
        query_len = len(self.orig_query)
        ref_len = len(self.orig_ref)
        print('\nQuery: {0} ({1} nt)'.format(self.q_name, query_len))
        print('Ref  : {0} ({1} nt)'.format(self.r_name, ref_len))
        print('\nQuery: {0: <7}{1} {2}'.format(self.q_pos + 1, self.aln[0][0], self.q_pos + 1 + query_len))
        print(' ' * 14 + aln_string)
        print ('Ref  : {0: <7}{1} {2}'.format(self.r_pos + 1, self.aln[0][1], self.r_pos + 1 + ref_len))
        print('\nScore: {0}\nCIGAR: {1}').format(self.score, MakeSamCigar(self.cigar))
        return

    def __str__(self):
        self.dump()
        return


class TestSuite(unittest.TestCase):

    # aln[0] is read
    # aln[1] is ref
    def test_matches(self):

        t = RealignReads()
        a1 = ('ACGTACGT', 'ACGTACGT')
        rd1 = 'ACGTACGT'  # perfect match
        self.assertEqual(t._MakeBamCigar(a1, rd1), [(0, 8)])
        a2 = ('ACTTTCGT', 'ACGTACGT')
        rd2 = 'ACTTTCGT'  # mismatches in middle
        self.assertEqual(t._MakeBamCigar(a2, rd2), [(0, 8)])
        a3 = ('CCGTACGT', 'ACGTACGT')
        rd3 = 'CCGTACGT'  # 1 mismatch at ref (subject) left end
        self.assertEqual(t._MakeBamCigar(a3, rd3), [(4, 1), (0, 7)])
        a4 = ('ACGTACGC', 'ACGTACGT')
        rd4 = 'ACGTACGC'  # 1 mismatch at read (query) left end
        self.assertEqual(t._MakeBamCigar(a4, rd4), [(0, 7), (4, 1)])
        a5 = ('CTGTACGT', 'ACGTACGT')
        rd5 = 'CTGTACGT'  # 2 mismatches at ref (subject) left end
        self.assertEqual(t._MakeBamCigar(a5, rd5), [(4, 2), (0, 6)])
        a6 = ('ACGTACCC', 'ACGTACGT')
        rd6 = 'ACGTACCC'  # 2 mismatches at read (query) left end
        self.assertEqual(t._MakeBamCigar(a6, rd6), [(0, 6), (4, 2)])

    def test_indels(self):

        t = RealignReads()
        a1 = ('ACGTACGT', 'ACGTACGT')
        rd1 = 'ACGTACGT'  # perfect match
        self.assertEqual(t._MakeBamCigar(a1, rd1), [(0, 8)])
        a2 = ('ACTTTCGT', 'ACGTACGT')
        rd2 = 'ACTTTCGT'  # mismatches in middle
        self.assertEqual(t._MakeBamCigar(a2, rd2), [(0, 8)])
        a3 = ('CCGTACGT', 'ACGTACGT')
        rd3 = 'CCGTACGT'  # 1 mismatch at ref (subject) left end
        self.assertEqual(t._MakeBamCigar(a3, rd3), [(4, 1), (0, 7)])
        a4 = ('ACGTACGC', 'ACGTACGT')
        rd4 = 'ACGTACGC'  # 1 mismatch at read (query) left end
        self.assertEqual(t._MakeBamCigar(a4, rd4), [(0, 7), (4, 1)])
        a5 = ('CTGTACGT', 'ACGTACGT')
        rd5 = 'CTGTACGT'  # 2 mismatches at ref (subject) left end
        self.assertEqual(t._MakeBamCigar(a5, rd5), [(4, 2), (0, 6)])
        a6 = ('ACGTACCC', 'ACGTACGT')
        rd6 = 'ACGTACCC'  # 2 mismatches at read (query) left end
        self.assertEqual(t._MakeBamCigar(a6, rd6), [(0, 6), (4, 2)])

    def test_softclips(self):

        t = RealignReads()
        a1 = ('ACGTACGT', '--GTACGT')
        rd1 = 'ACGTACGT'  # insertion to ref at the start
        self.assertEqual(t._MakeBamCigar(a1, rd1), [(4, 2), (0, 6)])
        a2 = ('ACGTACGT', 'ACGTAC--')
        rd2 = 'ACGTACGT'  # insertion to ref at the end
        self.assertEqual(t._MakeBamCigar(a2, rd2), [(0, 6), (4, 2)])
        a3 = ('ACGTACGT', '-CGTACG-')
        rd3 = 'ACGTACGT'  # insertion to ref at both ends
        self.assertEqual(t._MakeBamCigar(a3, rd3), [(4, 1), (0, 6), (4, 1)])
        a4 = ('ACGTAC--', 'ACGTACGT')
        rd4 = 'ACGTAC'  # deletion from the read at the end
        self.assertEqual(t._MakeBamCigar(a4, rd4), [(0, 6)])
        a5 = ('ACGTGG--', 'ACGTACGT')
        rd5 = 'ACGTGG'  # deletion from the read at the end
        self.assertEqual(t._MakeBamCigar(a5, rd5), [(0, 4), (4, 2)])
        '''
        a6 = ('--GTACGT', 'ACGTACGT')
        rd5 = 'ACGTACGT'  # check for a warning
        self.assertEqual(t._MakeBamCigar(a6, rd6), [(0, 6)])
        '''
        
    def test(self):
        import swalign as sw        
        ref = 'GTGATTAGTATTTTGTATTGTTGGGGTYGTTTTAGGGAGAAGATTTTTTTTTTTAT'
        rd1 = 'GTGATTAGTATTTGTATTGTTGGGGTTGTTTTAGGGAGAAGATTTTTTTNTTTTT'
        scoring_matrix = sw.ScoringMatrix('/home/jason/Dropbox/stopgap-measure/bismat.txt')
        
        sw = sw.LocalAlignment(scoring_matrix=scoring_matrix,
                               gap_penalty=-16,
                               gap_extension_penalty=-4,
                               gap_extension_decay=0.0,
                               prefer_gap_runs=True,
                               verbose=False,
                               globalalign=False)
        
        aln = sw.align(ref, rd1, ref_name='ref', query_name='rd1', rc=False)
        
        
        
     #def test_cigar_length(self):
     #    pass


if __name__ == '__main__':
    unittest.main()

'''
if __name__ == '__main__':
    s = ScoringMatrix('/home/ros259/Dropbox/stopgap-measure/bismat.txt')
    print(s.alphabet)
    print(s.scores)
    print(s.size)
    
'''

#matrix_file = '/home/jason/Dropbox/stopgap-measure/bismat.txt'
#n = realign.Needle(matrix=matrix_file)
'''
import nwalign as nw
nw.global_align('ACGT', 'ACGT', gap_open=-16,
                              gap_extend=-4,
                              matrix=matrix_file)
'''


