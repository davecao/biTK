# -*- coding: utf-8 -*-

import sys

__all__ = ['AlignInfo']

class AlignInfo(object):
    """
        Store aligned info
    """
    def __init__(self, query_name, target_name,
                    query_sequence, target_sequence, 
                    aligned_score, aligned_query, aligned_target,
                    method='Global', subMatrix='nuc42', tbMatrix=None,
                    gap_open=10.0, gap_extension=1.0,
                    gap_symbol='-', align_type='nt', t= None):
        """
            Alignment information for pairwised alignment

        Args:
            parameters (dict):
                gap_open (float) : 10.0
                gap_extension (float) : 1.0
                subMatrix (str) : 'nuc42' etc
                method (str): Global - Needleman-Wunsch
                align_type : 'nt' or 'pr'
                            nt - nucleotide
                            pr - protein
            query_name (str) : label/id of a query sequence
            target_name (str) : label/id of a target sequence
            query_sequence (str) : nucleotide/protein sequence
            target_sequence (str) : nucleotide/protein sequence

            n_match (int) : the number of matched pairs of query/target
            n_mismatch (int) : the number of mismatched pairs of query/target
            identities (float) : the percentage of matched pairs/non-gapped pairs
            positives  (float) : 
        """
        super(AlignInfo, self).__init__()
        # Input
        self.query_name = query_name
        self.target_name = target_name
        self.query_sequence = query_sequence
        self.target_sequence = target_sequence

        # Parameters for a pairwise alignment
        self.gap_open = gap_open
        self.gap_extension = gap_extension
        self.gap_symbol = gap_symbol
        self.subMatrix = subMatrix
        self.align_type = align_type
        self.method = method

        # Aligned results
        self.aligned_query = aligned_query
        self.aligned_target = aligned_target
        self.tbMatrix = tbMatrix

        # Statistics for an alignment
        self.aligned_score = aligned_score
        self.aligned_length = len(self.aligned_target)
        self.n_match = 0
        self.n_mismatch = 0
        self.identities = 0.0
        self.positives = 0.0
        self.n_gap_in_query = aligned_query.count(gap_symbol)
        self.n_gap_in_target = aligned_target.count(gap_symbol)
        # Support output formats
        self.fmt = ['txt', 'html', 'xml']
        self.o_fmt = 0

        # date:
        self.dateInfo = t
        # do statistics 
        self._do_stats()

    def _do_stats(self):
        """
            count match/mismatch
        """
        aligned_length = self.aligned_length
        for inx in range(aligned_length):
            t_item = self.aligned_target[inx]
            q_item = self.aligned_query[inx]
            if (t_item != self.gap_symbol) and (q_item != self.gap_symbol):
                # pairs is non gap
                if (t_item == q_item):
                    self.n_match += 1
                self.positives += 1
            elif (t_item == self.gap_symbol) or (q_item == self.gap_symbol):
                # gap vs nt/aa
                self.n_mismatch += 1
            else:
                # gap vs gap
                self.n_mismatch += 1
        self.identities = self.n_match/float(aligned_length) * 100
        self.positives = self.positives/aligned_length * 100

    def  _get_ofmt(self, fmt):
        """ 
            Return max value and its index in a list

        """
        #slist = [score_diag, score_up, score_left]
        #fmt_inx = [i for i, j in enumerate(self.fmt) if j==fmt]
        return self.fmt[fmt]

    def _output_txt(self):
        o_header ="""
###############################################
#Pairwised sequence alignment by {} 
# date : {}
# parameters:
#     gap_open: {}
#     gap_extension: {}
#     substitution matrix: {}
# Statistics:
#     raw score : {}
#     aligned length : {}
#     gaps in query: {}
#     gaps in target : {}
#     matched : {}
#     mismatched : {}
#     identities: {:.3f}
#     positives : {:.3f}
##############################################
""".format(self.method, self.dateInfo, self.gap_open, 
            self.gap_extension, 
            self.subMatrix, self.aligned_score, self.aligned_length,
            self.n_gap_in_query,self.n_gap_in_target,
            self.n_match, self.n_mismatch, 
            self.identities, self.positives)
        foldNumber=80
        print(o_header)
        print("{:10s}: {}".format(self.query_name, self.aligned_query))
        print("{:10s}: {}".format(self.target_name,self.aligned_target))

    def _output_html(self):
        pass
    def _output_xml(self):
        pass

    def print_tbmat(self):
        self.tbMatrix.print_dpMat()

    def show(self, fmt=0):
        """ 
            Show results of a pairwised sequence alignment 

        Args:

        Kwargs:
            format(int): 
                0 : text format 
                1 : html
                2 : xml
        """
        #fmt = self._get_ofmt(format.lower())
        #fmt_fun = eval('self._output_'+fmt)
        if fmt == 0:
            self._output_txt()
        elif fmt == 1:
            self._output_html()
        elif fmt == 2:
            self._output_xml()
        else:
            print("Error: Unknown specified output format {}".format(fmt))
            print("  fmt should be 0, 1 or 2.")
            sys.exit(1)

