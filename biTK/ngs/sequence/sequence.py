# -*- coding: utf-8 -*-
"""This module defines the sequence class, revised from weblogo library
"""
import os
import sys
import re
import math
import copy
import types
import multiprocessing.reduction

from biTK.ngs.sequence import nucleic_alphabet, phred33_alphabet, \
                            phred64_alphabet, Seq
from biTK.ngs.SubMatrice import SubstitutionMatrix as SUBMat
from biTK.ngs.sequence.align import Alignment 
from biTK import PY3K

if PY3K:
    string_types = str,
    integer_types = int,
    class_types = type,
    text_type = str
    binary_type = bytes
else:
    string_types = basestring,
    integer_types = (int, long)
    class_types = (type, types.ClassType)
    text_type = unicode
    binary_type = str

__all__ = ['Sequence', 'PHRED_ALPHABET']

PHRED_ALPHABET = {
    'phred33': phred33_alphabet,
    'phred64': phred64_alphabet
}
class _SequenceBuilderWrapper(object):
    def __call__(self, cls, state, **kwargs):
        print(state)
        inst = _SequenceBuilderWrapper()
        incls = cls()
        incls.__setstate__(state)
        inst.__class__ = incls
        return inst

class Sequence(object):
    """ 
    Sequence class for single protein sequence inherited from str class

    usage:
    >>> import biTK
    >>> from biTK.ngs.sequence import Alphabet, generic_alphabet, protein_alphabet, nucleic_alphabet
    >>> pr=biTK.ngs.sequence.Alphabet('ACSDGF')
    >>> pr='ACSDGF'
    >>> s1=biTK.ngs.sequence.Sequence(pr,name='test',alphabet=protein_alphabet)
    >>> print(repr(s1))
    >test
    ACSDGF
    """
    __slots__ = ['name', 'description', 'optionID','format', 'GC_content',
                'alphabet', 'quality_alphabet', 'quality_format',
                'raw_seq', 'quality_seq','scores', 'p_errors']

#    def __new__(cls, *args, **kwargs):
#        print("Sequence __new__ is executed")
#        #print(args)
#        #print(kwargs)
#        if not kwargs:
#            print("kwargs is empty.")
#        raw_seq, quality_seq = args
#        self = object.__new__(cls)
#        self.__init__(raw_seq, quality_seq, **kwargs)
#        return self

    def __init__(self, raw_seq, quality_seq, **kwargs):
    #def __init__(self, *args, **kwargs):
    #def __init__(self, raw_seq, quality_seq, alphabet, name, 
    #                   optionID, description, qs_fmt):
        """
        Store a sequence in fastq format
        
        Args:
            raw_seq(str)    : raw nt sequence
            quality_seq(str): quality sequence

        Kwargs: 
            alphabet (obj) : an object of Alphabet, default is generic_alphabet.
            name (str) : id/name of the sequence, default is None.  
            optionID (str) : optional Id of the sequence, default is None.
            description (str) : description of the sequence, default is None
            qs_fmt (str) : 'phred64' or 'phred33'.
            format (str) : 'sanger'- currently only used as label.
        """
        super(Sequence, self).__init__()
        print("{},{}".format(raw_seq, quality_seq))
        sys.stdout.flush()
        print(kwargs)
        sys.stdout.flush()

        #raw_seq = args[0]
        #quality_seq = args[1]
        qs_fmt = kwargs.get('qs_fmt', 'phred33')
        alphabet = kwargs.get('alphabet', nucleic_alphabet)
#        quality_alphabet (obj) : default is phred33_fmt.
        quality_alphabet = PHRED_ALPHABET[qs_fmt]

        self.alphabet = alphabet
        self.quality_format = qs_fmt
        self.quality_alphabet = quality_alphabet
        self.name = kwargs.get('name', None)
        # Check missing values
        if len(raw_seq) != len(quality_seq):
            print("Missing values:")
            print("The length of reads is not equal to that of quality scores")
            sys.exit(1)
        try:
            self.raw_seq = Seq(raw_seq, alphabet=alphabet)
            #print("{}:{}".format(quality_seq,quality_alphabet))
            #print("{} - {}:{}".format(self.name, raw_seq, quality_seq))
            #sys.stdout.flush()
        except ValueError:
            print("Failed to construct objects of class Seq")
            sys.exit(1)
        try:
            self.quality_seq = Seq(quality_seq, alphabet=quality_alphabet)
            #print("{}:{}".format(quality_seq,quality_alphabet))
            #print("{} - {}:{}".format(self.name, raw_seq, quality_seq))
            #sys.stdout.flush()
        except ValueError:
            print("Failed to construct objects of class Seq")
            sys.exit(1)
        #self.name = kwargs.get('name', None)
        self.description = kwargs.get('description', None)
        self.optionID = kwargs.get('optionID', None)
#        self.name = name
#        self.description = description
#        self.optionID = optionID
        self.scores, self.p_errors = self.__chr2score(format=qs_fmt)
        self.GC_content = self.__get_GC_content()
        

    @property
    def __deepcopy_keyattrs(self):
        return {'alphabet':self.alphabet, 
                'quality_alphabet':self.quality_alphabet,
                'name':self.name,
                'optionID':self.optionID,
                'quality_format':self.quality_format,
                'description':self.description
                }
    @property
    def __deepcopy_args(self):
        return [self.raw_seq.tostring(), self.quality_seq.tostring()]
    
    def __deepcopy__(self, memo):
        """ Copy object """
        kwds = self.__deepcopy_keyattrs
        args = self.__deepcopy_args
        cls = self.__class__
        result = cls.__new__(cls, *args, **kwds)
        memo[id(self)] = result
        result.__init__(*args, **kwds)
        #for k, v in self.__dict__.items():
        #    setattr(result, k, copy.deepcopy(v, memo))
        return result

    def __chr2score(self, format='phred33'):
        """
            Convert chr to score
        Kwargs:
            format (str): Quality format
                phred64 | phred33
        """
        p_errors = []
        scores = []
        cons = -0.1
        for score in self.quality_seq.ords():
            if score == 255:
                print("Specified score format seem not to be {}".format(format))
                print("{}\n{}".format(self.raw_seq, self.quality_seq))
                sys.exit(1)
            scores.append(score)
            p_errors.append(math.pow(10, score*cons))
        return (scores, p_errors)

    def __get_GC_content(self):
        """ Return GC content of a sequence """
        nt_G = 0.0
        nt_C = 0.0
        nt_all = dict(self.raw_seq.word_count(1, self.alphabet))
        s = sum(nt_all.values(), 0.0)
        if s == 0.0:
            s = 1.0
        if 'G' in nt_all:
            nt_G = nt_all['G']
        if 'C' in nt_all:
            nt_C = nt_all['C']
        gc = (nt_G + nt_C)/s*100
        return gc

    def __str__(self):
        return self.raw_seq

    def __repr__(self):
        return "{}".format(self.raw_seq)

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.raw_seq == other.raw_seq
        return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.__repr__())

    def _get_slots(self):
        all_slots = (getattr(cls, '__slots__', ()) for cls in self.__class__.__mro__)
        r = set(slot for slots in all_slots for slot in slots)
        return r

    def __getstate__(self):
        # subclass do not must be slots
        try:
            state = vars(self).copy()
        except TypeError:
            state = {}

        for slot in self._get_slots():
            try:
                val = getattr(self, slot)
                state[slot] = val
            except AttributeError:
                pass
        return state

    def __setstate__(self, state):
        for k in state:
            setattr(self, k , state[k])

    def __reduce__(self):
        #print("Sequence __reduce__")
        from copy_reg import __newobj__
#        if hasattr(self, '__getnewargs__'):
#            args = self.__getnewargs__()
#        else:
#            args = ()
#
#        if hasattr(self, '__getstate__'):
#            state = self.__getstate__()
#        elif hasattr(type(self), '__slots__'):
#            state = self.__dict__, {k: getattr(self, k) for k in type(self).__slots__}
#        else:
#            state = self.__dict__
#
        if isinstance(self, list):
            listitems = self
        else:
            listitems = None

        if isinstance(self, dict):
            dictitems = self.items()
        else:
            listitems = None
        args = (self.raw_seq.tostring(), self.quality_seq.tostring())
        state = self.__getstate__()

        return __newobj__, (type(self),)+args, state, None, None
        # for pickle
#        kwargs = {'alphabet':self.alphabet, 
#                'quality_alphabet':self.quality_alphabet,
#                'name':self.name,
#                'optionID':self.optionID,
#                'quality_format':self.quality_format,
#                'description':self.description
#                }
#        args = (self.raw_seq.tostring(), self.quality_seq.tostring())
#        cls = self.__class__
###        # without kwargs is OK
##        return cls, args, kwargs
#        #return (_SequenceBuilderWrapper(), (Sequence, object), state)
#        return (cls, (args[0], args[1],), kwargs,)
#        #return (cls, (self.raw_seq, self.quality_seq, self.alphabet, self.name, self.#optionID, self.description, self.quality_format, ))

    def __trim__(self, pos, start=0):
        """
            Trim n elements from the lists
        Args:
            pos (int) : the specified index of the list
            start (int) : 0 or 1
                0 - [pos:], trimming from the head 
                1 - [:pos], trimming from the tail
        """
        inx = { 0: {'s':pos,'e':None}, 1:{'s':None,'e':pos}}
        if start not in inx:
            print("Arguments Error: start should be 0 or 1.")
            sys.exit(1)
        self.raw_seq = self.raw_seq[inx['s']:inx['e']]
        self.quality_seq = self.quality_seq[inx['s']:inx['e']]
        self.scores = self.scores[inx['s']:inx['e']]
        self.p_errors = self.p_errors[inx['s']:inx['e']]
        self.GC_content = self.__get_GC_content()

    def __remove_adapter(self, adapter):
        """
            Remove adapter sequence from Sequence
            under the development 
        Args
            adapter (str): an adapter sequence
        """
        pass


    def toFastq(self):
        fastq_fmt = "@{}{}".format(self.description, os.linesep)
        fastq_fmt += "{}{}".format(self.raw_seq, os.linesep)
        fastq_fmt += "+{}{}".format(self.optionID, os.linesep)
        fastq_fmt += "{}{}".format(self.quality_seq, os.linesep)
        return fastq_fmt

    def tally(self, alphabet=None):
        return self.raw_seq.tally(alphabet=alphabet)

    def length(self):
        """ Return the length of raw sequence """
        return len(self.raw_seq)

    def get_raw_seq(self):
        """ Return raw sequence """
        return self.raw_seq

    def get_GC_content(self):
        """ Return GC content of a sequence """
#        nt_G = 0.0
#        nt_C = 0.0
#        nt_all = dict(self.raw_seq.word_count(1, alphabet))
#        s = sum(nt_all.values(), 0.0)
#        if s == 0.0:
#            s = 1.0
#        if 'G' in nt_all:
#            nt_G = nt_all['G']
#        if 'C' in nt_all:
#            nt_C = nt_all['C']
#        gc = (nt_G + nt_C)/s*100
        return self.GC_content

    def trim_score_head(self, score):
        """
            Trim n letters with lower score
        """
        #from bisect import bisect_right
        #inx = bisect_right(self.scores, score)
        success = False
        l = len(self.scores)
        #inx = next((self.scores.index(n) for n in self.scores if n > score), l)
        inx = next((n for n in range(l) if self.scores[n] > score), l)
        if inx <= l and inx != 0:
            # Trim inx bases from the head
            self.__trim__(inx, start=0)
            success = True
        else:
            # 1. out of range, i.e., all elements are less than or equal to score
            # should give some warns? or be quiet?
            # 2. No elements need to be trimmed.
            pass
        return success

    def trim_score_tail(self, score):
        """
            Trim n letters with lower score
        """
        success = False
        l = len(self.scores)
        inx = next((n for n in range(l-1, -1, -1) if self.scores[n] > score), l)
        if inx <= l:
            self.__trim__(inx, start=l)
            success = True
        else:
            # 1. out of range, i.e., all elements are less than or equal to score
            # should give some warns? or be quiet?
            # 2. No elements need to be trimmed.
            pass
        return success

    def trim_by_score(self, score, order=[0, 1]):
        """ Caution: there may be empty after trimming from one of ends
        Args:
            score (float): the threshold for triming elements in the list
            order (list) : [0, 1] or [1, 0]
                [0, 1] - trim bases from the head first, then the tail, default.
                [1, 0] - trim bases from the tail first, then the head
        """
        action = {
            0 : self.trim_score_head,
            1 : self.trim_score_tail
        }
        success = False
        if action[order[0]](score) and action[order[1]](score):
            success = True
        return success

    def trim_len_head(self, n):
        """
            Trim n bases from the head 
        """
        if n > len(self.raw_seq):
            return
        self.__trim__(n, start=0)

    def trim_len_tail(self, n):
        """
            Trim n bases from the tail
        """
        # Caution: python's list is zero-based, subscript 
        if n > len(self.raw_seq):
            return
        self.__trim__(n, start=1)

    def get_n_head(self, n):
        return self.raw_seq.first_n_el(n)

    def words(self, k, alphabet=None):
        return self.raw_seq.words(k, alphabet=alphabet)

    def word_count(self, k, alphabet=None):
        return self.raw_seq.word_count(k, alphabet=alphabet)

    def align(self, adapterSeq, algo='NWAlign',score_mat=None, align_type='nr',
                    gap_open=10.0, gap_ext=1.0, gap_symbol='-'):
        """
            Sequence alignment
        Args:
            adapterSeq, 
        Kwargs:
            algo (str) : 'NWAlign' - global pairwised alignment  
                         'SWAlign' - local pairwised alignment
            score_mat (str) : 'nuc42', 'nuc44' etc. 
                      see biTK.ngs.SubMatrice.SubstitutionMatrix
            gap_open (float) : default is 10.0.  
            gap_ext  (float) : default is 1.0.
        """
        if not isinstance(adapterSeq, Seq):
            print("Error: wrong input argument")
            sys.exit(1)

        # pairwised sequence alignment
        aligner = Alignment(algo='NWAlign', align_type="nr")
        alignInfo = aligner(self.raw_seq, adapterSeq, 
                            substitutionMatrix=score_mat,
                            gap_open=gap_open, gap_extension=gap_extension, 
                            gap_symbol=gap_symbol)
        query = alignInfo.aligned_query
        target = alignInfo.aligned_target 
        return query, target

    def trimAdapter(self, adapters):
        """
            Trim adapter sequence(s)
        Args:
            adapters(str or list, tuple of str):
                adapter sequence or a list/tuple of adapter sequences 
        """
        from itertools import chain
        if isinstance(adapters, string_types):
            self.__remove_adapter(adapters)
        elif isinstance(adapters, (list, tuple)):
            for adapter in chain(adapters):
                self.__remove_adapter(adapter)
